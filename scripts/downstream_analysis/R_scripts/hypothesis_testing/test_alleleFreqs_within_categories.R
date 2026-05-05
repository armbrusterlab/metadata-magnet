library(dplyr)
library(tidyr)
library(rstatix)

test_alleleFreqs_within_categories <- function(df, site, residue, group_var = "category", seq_colname = "sequence", adjust = "fdr") {
  # Made a more concise version.
  # df[[seq_colname]] is assumed to contain aligned sequences
  # site is assumed to be a single int, pre-adjusted upon entry into the GUI
  adj_site <- site

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

  # make a version of the table to print
  cat("=== ALLELE COUNTS ===\n")
  tab_withTotals <- as.table(rbind(observed, not_observed, observed+not_observed))
  dimnames(tab_withTotals) <- list(allele_observed = c("yes", "no", "Total"), group = groups)
  print(tab_withTotals)
  
  ### NEW: Calculate overall statistics for enrichment tests
  total_sequences <- sum(xtab)
  total_allele <- sum(xtab["yes", ])
  overall_prop <- total_allele / total_sequences
  
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("Total sequences:", total_sequences, "\n")
  cat("Total", residue, "alleles:", total_allele, "\n")
  cat(glue("Overall proportion: {round(overall_prop, 4)} ( {round(overall_prop * 100, 1)}% )"), "\n\n")  
  
  # Overall test - is there any significant variation?
  #cat("\nFisher's exact test:\n")
  p <- tryCatch(
    {
      test <- fisher_test(xtab, detailed = TRUE)
      test$p
    },
    error = function(e) {
        # Error handling
        #message("Error: ", e$message)
        cat("Dataset may be too large to run Fisher's exact test with exact p-value; running with simulate.p.value = TRUE\n")
        test <- fisher_test(xtab, detailed = TRUE, simulate.p.value = TRUE, B=2000) # for now, use the default value of B, but B=100000 seems reasonable
        test$p
    }
  )

  print(glue("Fisher's exact test: p={p}"))
  #print(glue("Fisher's exact test: p={format(p, digits = 20, scientific = FALSE)}"))
  
  if (p < 0.05) {
    cat("p from Fisher's exact test is significant (p < 0.05); conducting further comparisons.\n\n")
    # cat("Category-wise enrichment tests:\n") # print later
    enrichment_results <- data.frame(
      category = groups,
      total_sequences = colSums(xtab),
      observed_alleles = xtab["yes", ],
      expected_alleles = colSums(xtab) * overall_prop,
      observed_prop = xtab["yes", ] / colSums(xtab),
      p_value_fisher_2x2 = NA
    )
    
    # Perform tests for each category
    for(i in 1:nrow(enrichment_results)) {    
      # 2x2 Fisher's exact test for this category vs all others
      category_A <- enrichment_results$observed_alleles[i]
      category_nonA <- enrichment_results$total_sequences[i] - category_A
      other_A <- total_allele - category_A
      other_nonA <- (total_sequences - enrichment_results$total_sequences[i]) - other_A
      
      fisher_2x2 <- fisher.test(
        matrix(c(category_A, category_nonA, other_A, other_nonA), nrow = 2),
        alternative = "greater"
      )
      enrichment_results$p_value_fisher_2x2[i] <- fisher_2x2$p.value
    }
    
    # Adjust p-values for multiple testing
    enrichment_results$p_adj_fisher_2x2 <- p.adjust(enrichment_results$p_value_fisher_2x2, method = adjust)
    
    # Add significance flags
    enrichment_results$significant_fisher_2x2 <- enrichment_results$p_adj_fisher_2x2 < 0.05
    
    # Format and print results
    enrichment_results_formatted <- enrichment_results |>
      mutate(
        observed_prop = round(observed_prop, 4),
        expected_alleles = round(expected_alleles, 2),
        p_value_fisher_2x2 = format.pval(p_value_fisher_2x2, digits = 3),
        p_adj_fisher_2x2 = format.pval(p_adj_fisher_2x2, digits = 3)
      ) |>
      select(category, observed_alleles, expected_alleles, p_adj_fisher_2x2, significant_fisher_2x2) |>
      rename(p_adj = p_adj_fisher_2x2, significant = significant_fisher_2x2)
    
    row.names(enrichment_results_formatted) <- NULL
    # print(enrichment_results_formatted) # Print later

    # Pairwise tests - which specific pairs are different?
    # cat("\nPairwise Fisher's tests:") # print later
    pairwise_results <- pairwise_fisher_test(xtab, p.adjust.method = adjust, alternative = "greater")
    # print(pairwise_results) # print later

    # # Identify significantly different pairs
    # sig_pairs <- pairwise_results |> filter(p.adj < 0.05)
    # if(nrow(sig_pairs) > 0) {
    #   cat("\nShowing only significantly different category pairs (adjusted p < 0.05):\n")
    #   print(sig_pairs, n=nrow(sig_pairs))
    # } else {
    #   cat("\nNo significantly different category pairs after adjustment\n")
    # }

    ### Summary of findings
    cat("\n=== SUMMARY ===\n")
    cat("Categories significantly enriched for", residue, "allele:\n")
    enriched_categories <- enrichment_results |> 
      filter(significant_fisher_2x2)
    if(nrow(enriched_categories) > 0) {
      for(i in 1:nrow(enriched_categories)) {
        cat(glue("  -{enriched_categories$category[i]} (observed: {enriched_categories$observed_alleles[i]}/{enriched_categories$total_sequences[i]} = {round(enriched_categories$observed_prop[i] * 100, 1)}%; expected: {round(enriched_categories$expected_alleles[i], 1)}, p-value adjusted by {adjust} is (p = {formatC(enriched_categories$p_adj_fisher_2x2[i], format = 'e', digits = 2)}))\n\n"))
            
        # additional print statements to facilitate data entry
        # cat(glue("{enriched_categories$category[i]} (p = {formatC(enriched_categories$p_adj_fisher_2x2[i], format = 'e', digits = 2)}; {round(enriched_categories$observed_prop[i] * 100, 1)}%)"), "\n") 
        # cat(glue("expected {round(enriched_categories$expected_alleles[i] / enriched_categories$total_sequences[i] * 100, 1)}%"))
        # cat(glue("expected {round(overall_prop * 100, 1)}%"), "\n\n")
      }
    } else {
      cat("  None\n")
    }

    sig_pairs <- pairwise_results |> filter(p.adj < 0.05)
    if(nrow(sig_pairs) > 0) {
      cat("\nSignificantly different category pairs (adjusted p < 0.05):\n")
      print(sig_pairs, n=nrow(sig_pairs))
    } else {
      cat("\nNo significantly different category pairs after adjustment\n")
    }

    cat("\n=== FULL RESULTS ===\n")
    cat("Category-wise enrichment tests:\n")
    print(enrichment_results_formatted)

    cat("\nPairwise Fisher's tests:") 
    print(pairwise_results, n=nrow(pairwise_results))

  } else {
    cat("\nOverall p was not significant\n")
  }
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