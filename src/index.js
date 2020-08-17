import register from "higlass-register";

import GFFDataFetcher from "./GFFDataFetcher";

register(
  { dataFetcher: GFFDataFetcher, config: GFFDataFetcher.config },
  { pluginType: "dataFetcher" }
);
