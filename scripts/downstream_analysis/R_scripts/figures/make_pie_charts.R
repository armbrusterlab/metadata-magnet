library(tidyverse)
library(glue)
theme_set(theme_classic())


make_pie_charts <- function(df, outdir, group_var1, group_var2, width = NA, height = NA) {
  # inputs: metadata df (whatever the current df in the GUI is), names of two categorical columns (as strings)
  
  # prepare to label plots
  facet_counts <- df |>
    group_by(!!sym(group_var2)) |>
    summarize(n = n()) |>
    mutate(facet_label = paste0(.data[[group_var2]], "\n(n=", n, ")")) # facet_label is like group_var2 but with info added on group sizes

  df <- df |>
    left_join(facet_counts, by = group_var2)
    
  #print(head(df$facet_label))
  
  df |>
    group_by(!!sym(group_var2)) |>
    summarize(n = n()) |>
    print()

  p <- df |> 
    group_by(!!sym(group_var1), facet_label) |> summarize(n=n()) |>
    ggplot() + 
    geom_col(aes(x="", y=n, fill = !!sym(group_var1))) +
    coord_polar(theta = "y") +
    facet_wrap(vars(facet_label), scales="free") +
    labs(title = glue("Pie charts faceted by {group_var2}, colored by {group_var1}"), x="", y="")

  make_dir(outdir)
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

# creates directories with that name if they don't already exist
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}