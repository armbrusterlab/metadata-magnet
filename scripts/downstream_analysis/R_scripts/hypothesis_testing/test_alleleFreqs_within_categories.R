library(dplyr)
library(tidyr)
library(rstatix)

test_alleleFreqs_within_categories <- function(
    df, site, residue,
    group_var = "category",
    seq_colname = "sequence",
    adjust = "fdr",
    forceTest = "auto"   # "auto", "chisq", or "fisher"
) {
  # Validate forceTest
  if (!forceTest %in% c("auto", "chisq", "fisher")) {
    print(forceTest)
    stop('forceTest must be "auto", "chisq", or "fisher"')
  }
  
  adj_site <- site
  
  df |>
    group_by(!!sym(group_var)) |>
    summarize(count = n()) |>
    print()
  
  ### Build contingency table
  summary_df <- df |>
    mutate(allele = substr(df[[seq_colname]], adj_site, adj_site)) |>
    mutate(is_allele = (allele == residue)) |>
    group_by(!!sym(group_var), is_allele) |>
    summarize(n = n(), .groups = "drop") |>
    complete(!!sym(group_var), is_allele = c(TRUE, FALSE), fill = list(n = 0))
  
  groups <- summary_df |> filter(is_allele) |> pull(!!sym(group_var))
  observed <- summary_df |> filter(is_allele) |> pull(n)
  not_observed <- summary_df |> filter(!is_allele) |> pull(n)
  
  xtab <- as.table(rbind(observed, not_observed))
  dimnames(xtab) <- list(allele_observed = c("yes", "no"), group = groups)
  print(xtab)
  
  ### Overall stats
  total_sequences <- sum(xtab)
  total_allele <- sum(xtab["yes", ])
  overall_prop <- total_allele / total_sequences
  
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("Total sequences:", total_sequences, "\n")
  cat("Total", residue, "alleles:", total_allele, "\n")
  cat(glue("Overall proportion: {round(overall_prop, 4)} ({round(overall_prop*100,1)}%)\n"))
  
  ### Enrichment tests
  enrichment_results <- data.frame(
    category = groups,
    total_sequences = colSums(xtab),
    observed_alleles = xtab["yes", ],
    expected_alleles = colSums(xtab) * overall_prop,
    observed_prop = xtab["yes", ] / colSums(xtab),
    p_value = NA,
    test_used = NA
  )
  
  choose_test <- function(cont_table, forceTest) {
    expected <- suppressWarnings(chisq.test(cont_table, correct = FALSE)$expected)
    
    if (forceTest == "chisq") {
      return(list(name = "Chi-square",
                  p = chisq.test(cont_table, correct = TRUE)$p.value))
    }
    if (forceTest == "fisher") {
      return(list(name = "Fisher's exact",
                  p = fisher.test(cont_table, alternative = "greater")$p.value))
    }
    
    # Automatic mode
    if (all(expected >= 5)) {
      list(name = "Chi-square",
           p = chisq.test(cont_table, correct = TRUE)$p.value)
    } else {
      list(name = "Fisher's exact",
           p = fisher.test(cont_table, alternative = "greater")$p.value)
    }
  }
  
  ### Per-category tests
  for (i in 1:nrow(enrichment_results)) {
    category_A <- enrichment_results$observed_alleles[i]
    category_nonA <- enrichment_results$total_sequences[i] - category_A
    other_A <- total_allele - category_A
    other_nonA <- (total_sequences - enrichment_results$total_sequences[i]) - other_A
    
    cont_table <- matrix(c(category_A, category_nonA, other_A, other_nonA), nrow = 2)
    
    test <- choose_test(cont_table, forceTest)
    enrichment_results$test_used[i] <- test$name
    enrichment_results$p_value[i] <- test$p
  }
  
  enrichment_results$p_adj <- p.adjust(enrichment_results$p_value, method = adjust)
  enrichment_results$significant <- enrichment_results$p_adj < 0.05
  enrichment_results$direction <- NA # for chi-square tests
  
  for (i in 1:nrow(enrichment_results)) {
    if (enrichment_results$significant[i] &&
        enrichment_results$test_used[i] == "Chi-square") {
      
      obs <- enrichment_results$observed_alleles[i]
      exp <- enrichment_results$expected_alleles[i]
      
      if (obs > exp) {
        enrichment_results$direction[i] <- "enriched"
      } else if (obs < exp) {
        enrichment_results$direction[i] <- "lacking"
      } else {
        enrichment_results$direction[i] <- "no difference"
      }
    } else if (enrichment_results$significant[i] &&
            enrichment_results$test_used[i] == "Fisher") {
        enrichment_results$direction[i] <- "enriched"
    } else {
      enrichment_results$direction[i] <- "n/a"
    }
  }
  
  
  enrichment_results_formatted <- enrichment_results |>
    mutate(
      observed_prop = round(observed_prop, 4),
      expected_alleles = round(expected_alleles, 2),
      p_value = format.pval(p_value, digits = 3),
      p_adj = format.pval(p_adj, digits = 3)
    )
  
  row.names(enrichment_results_formatted) <- NULL
  print(enrichment_results_formatted)
  
  ### Overall association test
  full_expected <- suppressWarnings(chisq.test(xtab, correct = FALSE)$expected)
  
  cat("\n=== OVERALL ASSOCIATION TEST ===\n")
  
  overall_test <- choose_test(xtab, forceTest)
  cat("Using", overall_test$name, "\n")
  cat("Expected frequencies:\n")
  print(round(full_expected, 2))
  
  if (overall_test$name == "Chi-square") {
    test <- chisq.test(xtab, correct = FALSE)
    p <- test$p.value
    cat(glue("\nChi-square: X²={round(test$statistic,3)}, df={test$parameter}, p={format.pval(p,3)}\n"))
    
    if (p < 0.05) {
      cat("\nStandardized residuals:\n")
      print(round(test$stdres, 2))
    }
  } else {
    test <- fisher_test(xtab, detailed = TRUE)
    print(test)
    p <- test$p
  }
  
  ### Pairwise comparisons
  if (p < 0.05) {
    cat("\nOverall test significant; running pairwise comparisons\n")
    
    pairwise_results <- data.frame()
    n_groups <- length(groups)
    
    for (i in 1:(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        pair_table <- xtab[, c(i, j)]
        test <- choose_test(pair_table, forceTest)
        
        pairwise_results <- rbind(pairwise_results, data.frame(
          group1 = groups[i],
          group2 = groups[j],
          p = test$p,
          test_used = test$name
        ))
      }
    }
    
    pairwise_results$p_adj <- p.adjust(pairwise_results$p, method = adjust)
    pairwise_results$significant <- pairwise_results$p_adj < 0.05
    
    print(pairwise_results)
    
    sig_pairs <- pairwise_results |> filter(significant)
    if (nrow(sig_pairs) > 0) {
      cat("\nSignificant pairs:\n")
      print(sig_pairs)
    } else {
      cat("\nNo significant pairs after adjustment\n")
    }
  } else {
    cat("\nOverall test not significant\n")
  }
  
  ### Summary
  cat("\n=== SUMMARY ===\n")
  # cat(glue("Groups significantly enriched or lacking in {residue} allele:\n\n"))
  cat(glue("Groups significantly enriched in {residue} allele:\n\n"))
  
  enriched_categories <- enrichment_results |> filter(significant)
  
  # # could report on categories lacking the allele for chi-square test, but won't because Fisher test is one-sided as intended
  # if (nrow(enriched_categories) == 0) {
  #   cat("  None\n")
  # } else {
  #   for (i in 1:nrow(enriched_categories)) {
  #     dir <- enriched_categories$direction[i]
  #     
  #     cat("  -", enriched_categories$category[i],
  #         "(", dir, "; observed:", enriched_categories$observed_alleles[i], "/",
  #         enriched_categories$total_sequences[i], "=",
  #         round(enriched_categories$observed_prop[i] * 100, 1), "%; expected:",
  #         round(enriched_categories$expected_alleles[i], 1),
  #         ", test:", enriched_categories$test_used[i],
  #         ", p_adj =", formatC(enriched_categories$p_adj[i], format = "e", digits = 2),
  #         ")\n")
  #   }
  # }
  
  if (nrow(enriched_categories) == 0) {
    cat("  None\n")
  } else {
    for (i in 1:nrow(enriched_categories)) {
      dir <- enriched_categories$direction[i]
      
      if (dir == "enriched") {
        cat("  -", enriched_categories$category[i],
            "(observed:", enriched_categories$observed_alleles[i], "/",
            enriched_categories$total_sequences[i], "=",
            round(enriched_categories$observed_prop[i] * 100, 1), "%; expected:",
            round(enriched_categories$expected_alleles[i], 1),
            ", test:", enriched_categories$test_used[i],
            ", p_adj =", formatC(enriched_categories$p_adj[i], format = "e", digits = 2),
            ")\n")
      }
    }
  }
  
  invisible(list(
    contingency_table = xtab,
    enrichment_results = enrichment_results,
    overall_p_value = p,
    overall_test_used = overall_test$name,
    total_sequences = total_sequences,
    total_allele = total_allele,
    overall_proportion = overall_prop
  ))
}



calibrate_coords <- function(seq, sites_of_interest) {
  # seq is a reference sequence, and sites_of_interest are the sites in the reference seq you want to annotate in the seqlogo
  # in an alignment, gap characters are inserted, so the coordinates need to be calibrated accordingly
  # at this point assume sites_of_interest is a sorted numeric vector
  non_gap = 0
  i = 1 # R is 1-indexed
  last_site <- max(sites_of_interest)
  
  sites_adj <- c()
  
  while (non_gap < last_site & i <= nchar(seq)) {
    if (substr(seq, i, i) != "-") {
      non_gap = non_gap + 1
      
      # Only add to sites_adj when we're at a non-gap position AND it matches a site of interest
      if (non_gap %in% sites_of_interest) {
        sites_adj <- append(sites_adj, i)
      }
    } 
    
    i = i + 1
  }
  
  return(sites_adj)
}