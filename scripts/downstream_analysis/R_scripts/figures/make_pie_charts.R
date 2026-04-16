library(tidyverse)
library(glue)
theme_set(theme_classic())


make_pie_charts <- function(df, outdir, group_var1, group_var2, width = NA, height = NA) {
  # inputs: metadata df (whatever the current df in the GUI is), names of two categorical columns (as strings)
  p <- df |> 
    group_by(!!sym(group_var1), !!sym(group_var2)) |> summarize(n=n()) |>
    ggplot() + 
    geom_col(aes(x="", y=n, fill = !!sym(group_var1))) +
    coord_polar(theta = "y") +
    facet_wrap(vars(!!sym(group_var2)), scales="free") +
    labs(title = glue("Pie charts faceted by {group_var2}, colored by {group_var1}"), x="", y="")

  ggsave(
    filename = glue("{outdir}/pieChart_{group_var1}_{group_var2}.pdf"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
}

