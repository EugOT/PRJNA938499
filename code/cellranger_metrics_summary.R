library(readr)
metrics_summary_SRR23615077 <- read_csv("/data/PRJNA938499/cellranger/SRR23615077/outs/metrics_summary.csv") %>% mutate(Run = "SRR23615077")
metrics_summary_SRR23615078 <- read_csv("/data/PRJNA938499/cellranger/SRR23615078/outs/metrics_summary.csv") %>% mutate(Run = "SRR23615078")
metrics_summary_SRR23615081 <- read_csv("/data/PRJNA938499/cellranger/SRR23615081/outs/metrics_summary.csv") %>% mutate(Run = "SRR23615081")
metrics_summary_SRR23615082 <- read_csv("/data/PRJNA938499/cellranger/SRR23615082/outs/metrics_summary.csv") %>% mutate(Run = "SRR23615082")
metrics_summary <-
  bind_rows(
    metrics_summary_SRR23615077,
    metrics_summary_SRR23615078,
    metrics_summary_SRR23615081,
    metrics_summary_SRR23615082)

metrics_summary |>
  select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here("metrics_summary.tsv"))
