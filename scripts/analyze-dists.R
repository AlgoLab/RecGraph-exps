#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
fname <- "50genes-distperpair.5def-9all-gir-bwa.tsv"

if (length(args) == 1) {
  fname <- args[1]
}

print(paste("Analyzing", fname))

# install.packages("tidyverse")

library("readr")
library("dplyr")
library("tidyr")
library("knitr")

df <- read_tsv(fname, show_col_types = FALSE)

spec(df)
summary(df)

df.diff <- (
  df
  %>% pivot_longer(!c(Read, bwa), names_to = "tool", values_to = "ed")
  %>% mutate(difference = ed - bwa)
  %>% select(-c(`ed`, `bwa`))
  %>% pivot_wider(names_from="tool", values_from="difference")
)


# df.diff

## df.diff <- (
##   df
##   %>% mutate(
##         `rspoa-5.M2-X4-O4-E2` = `rspoa-5.M2-X4-O4-E2` - `bwa`,
##         `rspoa-9.R2-r0.1` = `rspoa-9.R2-r0.1` - `bwa`,
##         `rspoa-9.R4-r0.1` = `rspoa-9.R4-r0.1` - `bwa`,
##         `rspoa-9.R8-r0.1` = `rspoa-9.R8-r0.1` - `bwa`,
##         `giraffe` = `giraffe` - `bwa`,
##         )
##   %>% select(-c(`bwa`))
## )

df.diff.long <- (
  df.diff
  %>% mutate(abs.sum = abs(`rspoa-5.M2-X4-O4-E2`) +
               abs(`rspoa-9.R2-r0.1`) + abs(`rspoa-9.R4-r0.1`) + abs(`rspoa-9.R8-r0.1`) +
               abs(`giraffe`))
  %>% filter(abs.sum != 0)
  %>% filter(abs.sum <= 20)
  %>% select(-abs.sum)
  %>% pivot_longer(!Read, names_to = "tool", values_to = "difference")
  %>% select(-Read)
  %>% group_by(tool)
  %>% mutate(
        outlier.high = difference > quantile(difference, .75) + 1.50*IQR(difference),
        outlier.low = difference < quantile(difference, .25) - 1.50*IQR(difference)
      )
  %>% ungroup()
)


## g <- (
##   ggplot(
##     df.diff.long,
##     aes(x=tool, y=difference, fill=tool)
##   )
##   + geom_jitter(
##       color="black",
##       data = filter(df.diff.long, outlier.high ==T | outlier.low == T),
##       alpha= 0.75,
##       height=.05, width = .3)
##   + geom_boxplot(outlier.shape = NA, outlier.size=6)
##   + scale_y_continuous("Difference with bwa", breaks=-300:300, minor_breaks=NULL)
##   + theme_linedraw()
## )
## g

tools <- df.diff %>% select(-Read) %>% colnames() %>% sort()
(
  df.diff
  %>% pivot_longer(!Read, names_to = "tool", values_to = "difference")
  %>% select(tool, difference)
  %>% mutate(
        difference = case_when(
          difference > 5 ~ 1000,
          TRUE ~ difference,
        )
      )
  %>% count(difference, tool)
  %>% ungroup()
  %>% arrange(difference, tool)
  %>% pivot_wider(names_from = tool, values_from = n, values_fill = 0)
  %>% select(difference, all_of(tools))
  %>% knitr::kable()
)
## 1000 significa difference > 5


df.classif <- (
  df
  %>% mutate(
        classif.pp = case_when(
          `rspoa-5.M2-X4-O4-E2` == `bwa` ~ "b) OK",
          `rspoa-5.M2-X4-O4-E2` < `bwa` ~ "a) WRONG",
          `rspoa-5.M2-X4-O4-E2` > `bwa` ~ "c) SUBOPTIMAL",
          ),
        )
  %>% select(-giraffe)
  %>% mutate(pp=`rspoa-5.M2-X4-O4-E2`)
  %>% pivot_longer(!c(Read, classif.pp, bwa, pp), names_to = "tool", values_to = "ed")
  %>% mutate(
        classif.tool = case_when(
          `ed` == `bwa` ~ "b) OK",
          `ed` < `bwa` ~ "a) WRONG",
          `ed` > `bwa` ~ case_when(
                   `ed` < `pp` ~ "c1) SUBOPTIMAL (better)",
                   `ed` == `pp` ~ "c2) SUBOPTIMAL (==)",
                   `ed` > `pp` ~ "c3) SUBOPTIMAL (worse)",
                   )
          ),
      )
  %>% select(-c(`ed`, `bwa`, `pp`))
  )



tools <- df.classif$tool %>% unique() %>% sort()
(
  df.classif
  #%>% select(tool, difference)
  %>% count(classif.pp, classif.tool, tool)
  %>% ungroup()
  %>% pivot_wider(names_from = tool, values_from = n, values_fill = 0)
  %>% select(classif.pp, classif.tool, all_of(tools))
  %>% select(-`rspoa-5.M2-X4-O4-E2`)
  %>% knitr::kable()
)


(
  df
  #%>% select(Read, bwa, giraffe, `rspoa-5.M2-X4-O4-E2`)
  %>% pivot_longer(!c(Read, bwa), names_to = "tool", values_to = "ed")
  %>% mutate(
        classif.tool = case_when(
          `ed` == `bwa` ~ "b) OK",
          `ed` < `bwa` ~ "a) WRONG",
          `ed` > `bwa` ~ "c) SUBOPTIMAL",
          ),
      )
  %>% count(classif.tool, tool)
  %>% pivot_wider(names_from="tool", values_from="n")
  %>% knitr::kable()
)
