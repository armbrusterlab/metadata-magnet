library(dplyr)
library(tidyr)
library(rstatix)

test_alleleFreqs_within_categories <- function(df, site, residue, group_var = "category", seq_colname = "sequence", adjust = "fdr", forceChisq = FALSE) {
  # df[[seq_colname]] is assumed to contain aligned sequences
  # site is assumed to be a single int, pre-adjusted upon entry into the GUI
  adj_site <- site
  
  df |>
    group_by(!!sym(group_var)) |>
    summarize(count=n()) |>
    print()

  ### make contingency table for pairwise Fisher's test
  summary_df <- df |> 
    mutate(allele = substr(df[[seq_colname]], adj_site, adj_site)) |> 
    mutate(is_allele = (allele == residue)) |> 
    group_by(!!sym(group_var), is_allele) |> 
    summarize(n = n(), .groups = "drop") 
  
  # Ensure all combinations of group_var and is_allele exist
  summary_df <- summary_df |>
    complete(!!sym(group_var), is_allele = c(TRUE, FALSE), fill = list(n = 0))

  # get group names in the same order as in summary_df
  groups <- summary_df |> filter(is_allele == TRUE) |> pull(!!sym(group_var))
  observed <- summary_df |> filter(is_allele == TRUE) |> pull(n)
  not_observed <- summary_df |> filter(is_allele == FALSE) |> pull(n)
  
  xtab <- as.table(rbind(observed, not_observed))
  dimnames(xtab) <- list(allele_observed = c("yes", "no"), group = groups)
  print(xtab)
  
  ### Calculate overall statistics for enrichment tests
  total_sequences <- sum(xtab)
  total_allele <- sum(xtab["yes", ])
  overall_prop <- total_allele / total_sequences
  
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("Total sequences:", total_sequences, "\n")
  cat("Total", residue, "alleles:", total_allele, "\n")
  cat(glue("Overall proportion: {round(overall_prop, 4)} ( {round(overall_prop * 100, 1)}% )"), "\n")
  
  cat("\n=== ENRICHMENT TESTS ===\n")
  enrichment_results <- data.frame(
    category = groups,
    total_sequences = colSums(xtab),
    observed_alleles = xtab["yes", ],
    expected_alleles = colSums(xtab) * overall_prop,
    observed_prop = xtab["yes", ] / colSums(xtab),
    p_value = NA,
    test_used = NA
  )
  
  # Perform tests for each category
  for(i in 1:nrow(enrichment_results)) {    
    # 2x2 contingency table for this category vs all others
    category_A <- enrichment_results$observed_alleles[i]
    category_nonA <- enrichment_results$total_sequences[i] - category_A
    other_A <- total_allele - category_A
    other_nonA <- (total_sequences - enrichment_results$total_sequences[i]) - other_A
    
    cont_table <- matrix(c(category_A, category_nonA, other_A, other_nonA), nrow = 2)
    
    # Check expected frequencies for chi-square test
    expected <- chisq.test(cont_table, correct = FALSE)$expected
    # print("Expected frequencies for chi-square test:")
    # print(expected)
    
    # Choose test based on expected frequencies
    if (all(expected >= 5) | forceChisq) {
      # Use Pearson's chi-square test with continuity correction
      test_result <- chisq.test(cont_table, correct = TRUE)
      enrichment_results$test_used[i] <- "Chi-square"
      enrichment_results$p_value[i] <- test_result$p.value
    } else {
      # Fall back to Fisher's exact test
      test_result <- fisher.test(cont_table, alternative = "greater")
      enrichment_results$test_used[i] <- "Fisher's exact"
      enrichment_results$p_value[i] <- test_result$p.value
    }
  }
  
  # Adjust p-values for multiple testing
  enrichment_results$p_adj <- p.adjust(enrichment_results$p_value, method = adjust)
  
  # Add significance flags
  enrichment_results$significant <- enrichment_results$p_adj < 0.05
  
  # Format and print results
  enrichment_results_formatted <- enrichment_results |>
    mutate(
      observed_prop = round(observed_prop, 4),
      expected_alleles = round(expected_alleles, 2),
      p_value = format.pval(p_value, digits = 3),
      p_adj = format.pval(p_adj, digits = 3)
    )
  
  row.names(enrichment_results_formatted) <- NULL
  print(enrichment_results_formatted)
  
  # Determine which test to use for the overall test
  # Check expected frequencies for the full contingency table
  full_expected <- chisq.test(xtab, correct = FALSE)$expected
  
  cat("\n=== OVERALL ASSOCIATION TEST ===\n")
  if (all(full_expected >= 5) | forceChisq) {
    # Use Pearson's chi-square test
    cat("Using Pearson's Chi-square test\n")
    cat("Expected frequencies:\n")
    print(round(full_expected, 2))
    
    test <- chisq.test(xtab, correct = FALSE)
    cat("\nChi-square test results:\n")
    cat(glue("X-squared = {round(test$statistic, 3)}, df = {test$parameter}, p-value = {format.pval(test$p.value, digits = 3)}"))
    p <- test$p.value
    
    # Also show standardized residuals for significant results
    if (p < 0.05) {
      cat("\n\nStandardized residuals (|residual| > 2 indicate significant contribution):\n")
      residuals <- test$stdres
      print(round(residuals, 2))
    }
  } else {
    # Fall back to Fisher's exact test
    cat("Using Fisher's exact test (some expected frequencies < 5)\n")
    cat("Expected frequencies:\n")
    print(round(full_expected, 2))
    
    test <- fisher_test(xtab, detailed = TRUE)
    cat("\nFisher's exact test results:\n")
    print(test)
    p <- test$p
  }
  
  if (p < 0.05) {
    cat("\nOverall test is significant; conducting pairwise comparisons\n")
    
    # Pairwise tests - which specific pairs are different?
    cat("\n=== PAIRWISE COMPARISONS ===\n")
    cat(glue("Using {ifelse((all(full_expected >= 5) | forceChisq), 'chi-square', 'Fisher\\'s exact')} tests for pairwise comparisons\n"))
    
    # Get all pairwise combinations
    n_groups <- length(groups)
    pairwise_results <- data.frame()
    
    for(i in 1:(n_groups-1)) {
      for(j in (i+1):n_groups) {
        # Create 2x2 table for these two groups only
        pair_table <- xtab[, c(i, j)]
        
        # Check expected frequencies for this pair
        pair_expected <- chisq.test(pair_table, correct = FALSE)$expected
        
        # Choose test based on expected frequencies
        if (all(pair_expected >= 5) | forceChisq) {
          test_result <- chisq.test(pair_table, correct = TRUE)
          test_name <- "Chi-square"
          p_val <- test_result$p.value
        } else {
          test_result <- fisher.test(pair_table, alternative = "two.sided")
          test_name <- "Fisher's exact"
          p_val <- test_result$p.value
        }
        
        pairwise_results <- rbind(pairwise_results, data.frame(
          group1 = groups[i],
          group2 = groups[j],
          p = p_val,
          test_used = test_name,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Adjust p-values for multiple comparisons
    pairwise_results$p_adj <- p.adjust(pairwise_results$p, method = adjust)
    pairwise_results$significant <- pairwise_results$p_adj < 0.05
    
    # Print pairwise results
    print(pairwise_results)
    
    # Identify significantly different pairs
    sig_pairs <- pairwise_results |> filter(significant)
    if(nrow(sig_pairs) > 0) {
      cat("\nSignificantly different category pairs (adjusted p < 0.05):\n")
      print(sig_pairs)
    } else {
      cat("\nNo significantly different category pairs after adjustment\n")
    }
  } else {
    cat("\nOverall test was not significant\n")
  }
  
  ### Summary of findings
  cat("\n=== SUMMARY ===\n")
  cat(glue("Categories significantly enriched for {residue} allele (using {ifelse((all(full_expected >= 5) | forceChisq), 'Chi-square', 'Fisher\\'s exact')} enrichment tests):\n"))
  enriched_categories <- enrichment_results |> 
    filter(significant)
  if(nrow(enriched_categories) > 0) {
    for(i in 1:nrow(enriched_categories)) {
      cat("  -", enriched_categories$category[i], 
          "(observed:", enriched_categories$observed_alleles[i], "/", enriched_categories$total_sequences[i],
          "=", round(enriched_categories$observed_prop[i] * 100, 1), "%; expected:", 
          round(enriched_categories$expected_alleles[i], 1), 
          ", test:", enriched_categories$test_used[i],
          ", p-value adjusted by", adjust, "is (p =", formatC(enriched_categories$p_adj[i], format = "e", digits = 2), ") )\n")
          
      # additional print statements to facilitate data entry
      cat(glue("{enriched_categories$category[i]} (p = {formatC(enriched_categories$p_adj[i], format = 'e', digits = 2)}; {round(enriched_categories$observed_prop[i] * 100, 1)}%)"), "\n") 
      cat(glue("expected {round(overall_prop * 100, 1)}%"), "\n\n")
    }
  } else {
    cat("  None\n")
  }
  
  # Return all results for potential further analysis
  invisible(list(
    contingency_table = xtab,
    enrichment_results = enrichment_results,
    overall_p_value = p,
    overall_test_used = ifelse((all(full_expected >= 5) | forceChisq), "Chi-square", "Fisher's exact"),
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