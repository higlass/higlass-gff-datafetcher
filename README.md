# GFF Data Fetcher for HiGlass

Quickly load data from a GFF file to create a gene annotations track in HiGlass.

[![HiGlass](https://img.shields.io/badge/higlass-🌸-brightgreen.svg)](http://higlass.io)

This is the source code for multivec tracks in HiGlass only; for the rest of HiGlass,
see these repositories:

- HiGlass viewer: https://github.com/hms-dbmi/higlass
- HiGlass server: https://github.com/hms-dbmi/higlass-server
- HiGlass docker: https://github.com/hms-dbmi/higlass-docker

## Usage

The live scripts can be found at:

- https://unpkg.com/higlass-gff-datafetcher/dist/higlass-gff-datafetcher.min.js

Configure the track in your view config; you should be all set from here!

```
[...
  {
    "type": "horizontal-gene-annotations",
    "height": 80,
    "data": {
      "type": "gff",
      "url": "https://pkerp.s3.amazonaws.com/public/GCF_001461035.1_ASM146103v1_genomic.gff.gz",
      "chromSizesUrl": "https://domain.com/my.chrom.sizes",
      "options": {
        "namePaths": [
          "gene",
          "annotationName"
        ]
      }
    }
  }
]
```

Note that the `chromSizesUrl` option is optional and only needs to be provided if the gff file lacks the `region` entries listing the chromosomes in the assembly.

For an example, see [`src/index.html`](src/index.html).

## Development

### Testing

To run the test suite:

```
npm run test-watch
```

### Installation

```bash
$ git clone https://github.com/higlass/higlass-gff-datafetcher
$ cd higlass-gff-datafetcher
$ npm install
```

If you have a local copy of higlass, you can then run this command in the higlass-gff-datafetcher directory:

```bash
npm link higlass
```

### Commands

- **Developmental server**: `npm start`
- **Production build**: `npm run build`
