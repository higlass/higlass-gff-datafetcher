import slugid from "slugid";
import pako from "pako";
import gff from "@gmod/gff";
import jp from "jsonpath";
import { tsvParseRows } from 'd3-dsv';
import { text } from 'd3-request';

/**
 * Shuffles array in place.
 * @param {Array} a items An array containing the items.
 */
function shuffle(a) {
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    const x = a[i];
    a[i] = a[j];
    a[j] = x;
  }
  return a;
}

/**
 * Take a list of genes, which can be any list with elements containing
 * { start, end } fields and return another list of { start, end }
 * fields containing the collapsed genes.
 *
 * The segments should be sorted by their start coordinate.
 *
 * The scale parameter is the number of base pairs per pixels
 */
function collapse(segments, scale) {
  const collapsed = [];

  // the maximum distance we allow between segments before collapsing them
  const MAX_DIST_BETWEEN = 5;

  // no segments in, no segments out
  if (!segments.length) {
    return [];
  }

  // start with the first segment
  let currStart = segments[0].xStart;
  let currEnd = segments[0].xEnd;

  // continue on to the next segments
  for (let i = 1; i < segments.length; i++) {
    if (segments[i].xStart < currEnd + (MAX_DIST_BETWEEN * 1) / scale) {
      // this segment is within merging distance -- merge it
      currEnd = Math.max(currEnd, segments[i].xEnd);
    } else {
      // this segment is outside of the merging distance, dump the current
      // collapsed segment and start a new one
      collapsed.push({
        type: "filler",
        xStart: currStart,
        xEnd: currEnd,
      });

      // start a new collapsed segment
      currStart = segments[i].xStart;
      currEnd = segments[i].xEnd;
    }
  }

  // add the final segment
  collapsed.push({
    xStart: currStart,
    xEnd: currEnd,
    strand: segments[0].strand,
    fields: [],
    type: "filler",
    uid: slugid.nice(),
  });

  return collapsed;
}

const chrToAbs = (chrom, chromPos, chromInfo) =>{
  return chromInfo.chrPositions[chrom].pos + chromPos;
}

function parseChromsizesRows(data) {
  const cumValues = [];
  const chromLengths = {};
  const chrPositions = {};

  let totalLength = 0;

  for (let i = 0; i < data.length; i++) {
    const length = Number(data[i][1]);
    totalLength += length;

    const newValue = {
      id: i,
      chr: data[i][0],
      pos: totalLength - length,
    };

    cumValues.push(newValue);
    chrPositions[newValue.chr] = newValue;
    chromLengths[data[i][0]] = length;
  }

  return {
    chrToAbs: ([chrName, chrPos]) =>
      chrToAbs(chrName, chrPos, { chrPositions }),
    cumPositions: cumValues,
    chrPositions,
    totalLength,
    chromLengths,
  };
}

function ChromosomeInfo(filepath, success) {
  const ret = {};

  ret.absToChr = (absPos) => (ret.chrPositions ? absToChr(absPos, ret) : null);

  ret.chrToAbs = ([chrName, chrPos] = []) =>
    ret.chrPositions ? chrToAbs(chrName, chrPos, ret) : null;

  return text(filepath, (error, chrInfoText) => {
    if (error) {
      // console.warn('Chromosome info not found at:', filepath);
      if (success) success(null);
    } else {
      const data = tsvParseRows(chrInfoText);
      const chromInfo = parseChromsizesRows(data);

      Object.keys(chromInfo).forEach((key) => {
        ret[key] = chromInfo[key];
      });
      if (success) success(ret);
    }
  });
}

const GFFDataFetcher = function GFFDataFetcher(HGC, ...args) {
  if (!new.target) {
    throw new Error(
      'Uncaught TypeError: Class constructor cannot be invoked without "new"'
    );
  }

  function gffObjToChromsizes(gffObj) {
    const annotations = gffObj
      .filter((x) => x[0].type === "region")
      .map((x) => x[0]);

    const chromSizes = parseChromsizesRows(
      annotations.map((x) => [x.seq_id, x.end])
    );

    chromSizes.chrToAbs = ([chrom, pos]) => chrToAbs(chrom, pos, chromSizes);
    return chromSizes;
  }

  function gffToHgGene(gb, namePaths, chromSizes) {
    const importance = gb.end - gb.start;
    const uid = slugid.nice();

    let queryName = null;

    if (namePaths) {
      for (const namePath of namePaths) {
        queryName = jp.query(gb.attributes, namePath)[0];
        if (queryName) break;
      }
    }
    return {
      xStart: chromSizes.chrToAbs([gb.seq_id, gb.start]),
      xEnd: chromSizes.chrToAbs([gb.seq_id, gb.end]),
      strand: gb.strand,
      chrOffset: chromSizes.chrPositions[gb.seq_id].pos,
      importance: gb.end - gb.start,
      uid,
      type: gb.type,
      fields: [
        gb.seq_id,
        gb.start,
        gb.end,
        (queryName && queryName[0]) || "",
        importance,
        gb.strand,
        "",
        "",
        gb.type,
        gb.name,
        gb.start.toString(),
        gb.end.toString(),
        gb.start.toString(),
        gb.end.toString(),
      ],
    };
  }

  class GFFDataFetcherClass {
    constructor(dataConfig) {
      this.dataConfig = dataConfig;
      this.trackUid = slugid.nice();

      this.errorTxt = "";
      this.dataPromise = this.loadGff(dataConfig);
    }

    async loadGff(dataConfig) {
      if (dataConfig.chromSizesUrl) {
        await new Promise((resolve) => {
          ChromosomeInfo(dataConfig.chromSizesUrl, (chromInfo) => {
            this.chromSizes = chromInfo;
            resolve();
          });
        });
      }

      if (dataConfig.url) {
        const extension = dataConfig.url.slice(dataConfig.url.length - 3);
        const gzipped = extension === ".gz";

        try {

        const response = await fetch(dataConfig.url, {
          mode: "cors",
          redirect: "follow",
          method: "GET",
        });

        const buffer = gzipped ? await response.arrayBuffer()  : await response.text();
        const gffText = gzipped
          ? pako.inflate(buffer, { to: "string" })
          : buffer;
        // store all the GFF file annotations
        this.gffObj = gff.parseStringSync(gffText);

        this.createGenesAndChroms();
        } catch (err)  {
          console.error("err:", err);
        }
      } else if (dataConfig.text) {
        await new Promise((resolve, reject) => {
          this.gffObj = gff.parseStringSync(dataConfig.text);

          this.createGenesAndChroms();
          resolve();
        });
      } else {
        console.error(
          'Please enter a "URL" or a "text" field to the data config'
        );
      }
    }

    createGenesAndChroms() {
      const excludeTypes = this.dataConfig.options && this.dataConfig.options.excludeTypes;

      let items = this.gffObj.map(x => x[0]);

      if (excludeTypes) {
        items = items.filter(x => !excludeTypes.includes(x.type))
      }

      this.genes = items;

      if (!this.chromSizes) {
        this.chromSizes = gffObjToChromsizes(this.gffObj);
      }

      this.hgGenes = this.genes.map((x) =>
        gffToHgGene(
          x,
          this.dataConfig.options && this.dataConfig.options.namePaths,
          this.chromSizes
        )
      );
    }

    tilesetInfo(callback) {
      this.tilesetInfoLoading = true;

      return this.dataPromise
        .then(() => {
          this.tilesetInfoLoading = false;

          const TILE_SIZE = 1024;
          let retVal = {};

          const totalLength = this.chromSizes.totalLength;

          // retVal[this.trackUid] = {
          retVal = {
            tile_size: TILE_SIZE,
            max_zoom: Math.ceil(
              Math.log(totalLength / TILE_SIZE) / Math.log(2)
            ),
            max_width: totalLength,
            min_pos: [0],
            max_pos: [totalLength],
          };

          if (callback) {
            callback(retVal);
          }

          return retVal;
        })
        .catch((err) => {
          this.tilesetInfoLoading = false;

          console.error(err);

          if (callback) {
            callback({
              error: `Error parsing gff: ${err}`,
            });
          }
        });
    }

    fetchTilesDebounced(receivedTiles, tileIds) {
      const tiles = {};

      const validTileIds = [];
      const tilePromises = [];

      for (const tileId of tileIds) {
        const parts = tileId.split(".");
        const z = parseInt(parts[0], 10);
        const x = parseInt(parts[1], 10);

        if (Number.isNaN(x) || Number.isNaN(z)) {
          console.warn("Invalid tile zoom or position:", z, x);
          continue;
        }

        validTileIds.push(tileId);
        tilePromises.push(this.tile(z, x));
      }

      Promise.all(tilePromises).then((values) => {
        for (let i = 0; i < values.length; i++) {
          const validTileId = validTileIds[i];
          tiles[validTileId] = values[i];
          tiles[validTileId].tilePositionId = validTileId;
        }

        receivedTiles(tiles);
      });
      // tiles = tileResponseToData(tiles, null, tileIds);
      return tiles;
    }

    tile(z, x) {
      return this.tilesetInfo().then((tsInfo) => {
        const tileWidth = +tsInfo.max_width / 2 ** +z;

        // get the bounds of the tile
        const minX = tsInfo.min_pos[0] + x * tileWidth;
        const maxX = tsInfo.min_pos[0] + (x + 1) * tileWidth;

        const scaleFactor = 1024 / 2 ** (tsInfo.max_zoom - z);

        const filtered = this.hgGenes.filter(
          (v) => v.xEnd > minX && v.xStart < maxX
        );

        const collapsedPlus = collapse(
          filtered.filter((v) => v.strand === "+"),
          scaleFactor
        );
        const collapsedMinus = collapse(
          filtered.filter((v) => v.strand !== "+"),
          scaleFactor
        );

        collapsedPlus.forEach((v) => {
          v.strand = "+";
        });
        collapsedMinus.forEach((v) => {
          v.strand = "-";
        });

        const shuffledFiltered = shuffle(filtered);

        let values = [];
        const TILE_CAPACITY = 20;
        // fill the tile with entries that are within it
        for (let i = 0; i < shuffledFiltered.length; i++) {
          if (values.length >= TILE_CAPACITY) break;

          values.push(shuffledFiltered[i]);
        }

        values = [...values, ...collapsedPlus, ...collapsedMinus];
        // values = values.concat(collapsedPlus).concat(collapsedMinus);
        // we're not going to take into account importance

        return values;
      });
    }
  }

  return new GFFDataFetcherClass(...args);
}; // end function wrapper

GFFDataFetcher.config = {
  type: "gff",
};

export default GFFDataFetcher;
