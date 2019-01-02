gcloud \
  alpha genomics pipelines run \
  --pipeline-file wdl_pipeline.yaml \
  --inputs-from-file WDL=test-wdl/ga4ghMd5.wdl,\
WORKFLOW_INPUTS=test-wdl/ga4ghMd5.inputs.json,\
WORKFLOW_OPTIONS=test-wdl/basic.papi.us.options.json \
  --env-vars WORKSPACE=gs://vumc/biostat/cqs/test/workspace,\
OUTPUTS=gs://vumc/biostat/cqs/test/outputs \
  --logging gs://vumc/biostat/cqs/test/loggings/md5.txt
