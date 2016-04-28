library("plyr")
library("dplyr")
library("tidyr")

bin_up <- function(x, tdat2, unit = "Subject") {
    agg <- x %>%
        inner_join(tdat2, "TrialID") %>%
        group_by_(unit, quote(Group), quote(CompType),
                  quote(Cond), quote(FID), quote(ID)) %>%
        summarize(Y = n()) %>%
        ungroup()

    agg_n <- agg %>%
        group_by_(unit, quote(Group), quote(CompType),
                  quote(Cond), quote(FID)) %>%
        summarize(N = sum(Y)) %>%
        ungroup() %>%
        inner_join(agg,
                   c(unit, "Group", "CompType", "Cond", "FID"))

    all_frames <- agg_n %>%
        select_(unit, quote(Group), quote(CompType), quote(Cond), quote(FID)) %>%
        distinct() %>%
        merge(agg_n %>% select(ID) %>% distinct())

    all_frames %>%
        left_join(agg_n,
                  c(unit, "Group", "CompType", "Cond", "FID", "ID")) %>%
        as_data_frame() %>%
        mutate(Y = ifelse(is.na(Y), 0, Y),
               N = ifelse(is.na(N), 0, N),
               CompType = factor(CompType),
               Cond = factor(Cond))    
}

agg_over_units <- function(x, full = FALSE) {
    ff <- x %>%
        group_by(Group, CompType, Cond, FID, ID) %>%
        summarize(Y = sum(Y, na.rm = TRUE), N = sum(N, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(p = Y / N)
    if (full) ff else ff %>% `[[`("p")
}

boot_by_unit <- function(x, unit = "Subject") {
    ff <- split(x, f = x[[unit]])
    sample(ff, length(ff), replace = TRUE) %>% bind_rows()
}


agg_up2 <- function(x, full = FALSE) {
    ff <- x %>%
        group_by(Group, CompType, Cond, FID) %>%
        summarize(C = sum(C), N = sum(NC)) %>%
        ungroup() %>%
        mutate(p = C / (C + N))
    if (full) ff else ff %>% `[[`("p")
}

add_y <- function(rx, dat) {
    cat(unlist(rx), "\n")
    data_frame(y = dat %>%
                   filter(Group2 == rx[["Group2"]],
                          CompType == rx[["CompType"]],
                          between(ms, rx[["xmin"]], rx[["xmax"]])) %>%
                              `[[`("UL") %>% max())
}

fit_mnom <- function(x, form = cbind(C, NC) ~ G * CT * CM) {
    junk <- capture.output(
        mod <- nnet::multinom(form, x)
        )
    v <- as.vector(coef(mod), "numeric")
    names(v) <- colnames(coef(mod))
    v
}

permute_unit <- function(x, unit_codes, vname_short) {
    all_but_var <- x[setdiff(tbl_vars(x), vname_short)]
    join_fields <- intersect(tbl_vars(all_but_var),
                             tbl_vars(unit_codes))
    unit_codes[[vname_short]] <- sample(unit_codes[[vname_short]])
    unit_codes %>%
        inner_join(all_but_var, join_fields)
}

## ulevels are the unit levels
sync_perm <- function(x, cond_codes, var, ulevels) {
    all_units <- do.call("select_", c(list(x), as.list(ulevels))) %>%
        distinct()

    dat_split <- split(all_units,
                       lapply(ulevels[-1], function(vx) all_units[[vx]]))

    n_per <- dat_split %>% lapply(nrow) %>% unlist()

    min_unit <- n_per[which(n_per == min(n_per))[1]]

    svec1 <- sample(c(-1, 1), min_unit, TRUE)

    all_units[["svec"]] <-
        sapply(seq_along(n_per), function(ix) {
                   sample(c(rep(1, n_per[ix] - min_unit), svec1))
               }) %>% unlist()

    all_var <- inner_join(cond_codes, all_units, ulevels)


    all_var[[var]] <- all_var[[var]] * all_var[["svec"]]

    all_but_var <- x[setdiff(tbl_vars(x), var)]
    join_fields <- intersect(tbl_vars(all_var), tbl_vars(all_but_var))

    all_var %>%
        inner_join(all_but_var, join_fields)
}

bin_comp <- function(x, tdat2, unit) {
    ff <- x %>% inner_join(tdat2, "TrialID") %>%
        mutate(Comp = ifelse(ID == "Comp", "C", "NC")) %>%
        group_by_(unit, quote(Group), quote(CompType), quote(Cond),
                  quote(FID), quote(Comp)) %>%
        summarize(Y = n()) %>% ungroup()

    all_comp <- ff %>% select(Comp) %>% distinct()

    all_fr <- ff %>% select_(unit, quote(Group), quote(CompType),
                             quote(Cond), quote(FID)) %>% distinct() %>%
                                 merge(all_comp) %>% as_data_frame()

    join_fr <- intersect(tbl_vars(ff), tbl_vars(all_fr))
    all_fr %>%
        left_join(ff, join_fr) %>%
        mutate(Y = ifelse(is.na(Y), 0, Y)) %>%
        spread(Comp, Y)
}

sync_boot <- function(x, unit, iv) {
    udat <- x %>% select_(iv, unit) %>% distinct()
    dat_split <- split(udat, udat[iv])
    n_per <- dat_split %>% lapply(nrow) %>% unlist()
    min_n <- n_per[which(n_per == min(n_per))[1]]

    boot_ix <- sample(seq_len(min_n), min_n, TRUE)
    boot_n <- table(boot_ix) %>% as.integer()

    bdat <- lapply(dat_split, function(dx) {
               extra <- sample(seq_len(nrow(dx)), nrow(dx) - sum(boot_n),
                               replace = FALSE)
               slix <- c(extra, rep(sample(setdiff(seq_len(nrow(dx)), extra),
                          length(boot_n), FALSE), boot_n))
               slice(dx, slix)
           }) %>% bind_rows()

    bdat %>%
        inner_join(x, c(unit, iv))
}

get_boot_pvals <- function(x, eff_keep, unit, iv, nmc = 1000,
                           mod_form = cbind(C, NC) ~ G * CT * CM) {
    orig <- daply(x, .(FID), fit_mnom,
                  form = mod_form)[, eff_keep, drop = FALSE]

    ax <- replicate(nmc, x %>%
                        sync_boot(unit, iv) %>%
                        daply(.(FID), fit_mnom, form = mod_form))

    boot_sd <- apply(ax[, eff_keep, , drop = FALSE], c(1, 2), sd)

    t_val <- abs(orig / boot_sd)
    2 * (1 - pnorm(t_val)) * sign(orig)
}

get_boot_pvals_ws <- function(x, eff_keep, unit, nmc = 1000,
                           mod_form = cbind(C, NC) ~ G * CT * CM) {
    orig <- daply(x, .(FID), fit_mnom,
                  form = mod_form)[, eff_keep, drop = FALSE]

    ax <- replicate(nmc,
                    x %>%
                        boot_by_unit(unit) %>%
                        daply(.(FID), fit_mnom, form = mod_form)
                    )

    boot_sd <- apply(ax[, eff_keep, , drop = FALSE], c(1, 2), sd)

    t_val <- abs(orig / boot_sd)
    2 * (1 - pnorm(t_val)) * sign(orig)
}

get_clusters <- function(x, max_cms_only = FALSE) {
    tvec <- (abs(x) < .05) * sign(x)
    tvec_rle <- rle(tvec)
    run_sig <- tvec_rle$values != 0
    res <- data_frame()
    if (sum(run_sig) == 0) {
        if (max_cms_only) {res <- 0} else {}
    } else {
        t0 <- sapply(which(run_sig), function(cx) {
                         sum(tvec_rle$lengths[seq_len(cx - 1)]) + 1
                     })
        t1 <- mapply(function(ix, iy) {tvec_rle$lengths[ix] + iy - 1},
                     which(run_sig), t0, SIMPLIFY = FALSE) %>% unlist()
        cms <- mapply(function(ix, iy) {
                          sum(-2 * log(abs(x[ix:iy])))
                      }, t0, t1)
        names(t0) <- NULL
        clust <- data_frame(run_id = seq_len(sum(run_sig)),
                            t0 = names(x)[t0],
                            t1 = names(x)[t1],
                            cms = cms)
        res <- if (max_cms_only) max(cms) else clust
    }
    return(res)
}

get_p_value <- function(rx, pmx) {
    eff <- as.character(rx[["Effect"]][1])
    sapply(rx[["cms"]],
           function(x) sum(c(x, pmx[, eff]) >= x)) / (length(pmx[, eff]) + 1)
}

cluster_pvalues <- function(orig, pvals, pmx_1, pmx_2,
                            efflist = c("G", "CM", "CT:CM", "G:CT",
                                "G:CM", "G:CT:CM")) {
    calc_mean_paramest <- function(x, ori) {
        from_pm <- which(rownames(ori) == x[["t0"]])
        to_pm <- which(rownames(ori) == x[["t1"]])
        data_frame(mpe = mean(ori[, x[["Effect"]]][from_pm:to_pm]))
    }
    ori <- readRDS(orig)
    pval <- readRDS(pvals)
    pmx1 <- readRDS(pmx_1)

    pmx_full <- pmx1
    if (!is.null(pmx_2)) {
        pmx2 <- readRDS(pmx_2)
        pmx_full <- cbind(pmx1, pmx2)
    } else {}
    cols_keep <- intersect(intersect(colnames(ori), colnames(pval)), colnames(pmx_full))

    clust <- adply(pval[, cols_keep, drop = FALSE], 2, get_clusters) %>%
        rename(Effect = X1) %>% mutate(Effect = as.character(Effect))
    clust2 <- clust %>%
        group_by(Effect, run_id) %>% do(calc_mean_paramest(., ori)) %>%
        inner_join(clust, c("Effect", "run_id"))

    clust2 %>%
        filter(Effect %in% efflist) %>%
        group_by(Effect, run_id) %>%
        do(pval = get_p_value(., pmx_full)) %>% unnest() %>%
        inner_join(clust2)
}


## add 200 ms to account for EM delay (for adults)
tdat2 <- readRDS("derived/trial_data.rds") %>%
    as_data_frame() %>%
    mutate(frEnd = round(60 * ((DPlag + 200) / 1000)))

edat_ctt <- readRDS("derived/eye_data_cumulative.rds") %>%
    as_data_frame()

dp_frames <- edat_ctt %>%
    inner_join(tdat2, "TrialID") %>%
    rename(orig_FID = FID) %>%
    mutate(FID = orig_FID - frEnd) %>%
    filter(FID >= -31, FID <= 91) %>%
    select(TrialID, FID, orig_FID, ID, Pad)

dp_subj <- dp_frames %>%
    inner_join(tdat2, "TrialID") %>%
    mutate(Comp = ifelse(ID == "Comp", "C", "NC"),
           bin = floor((FID + 1) / 3) * 3) %>%
    count(Subject, Group, CompType, Cond, FID = bin, Comp) %>%
    ungroup() %>%
    spread(Comp, n) %>%
    mutate(C = ifelse(is.na(C), 0, C),
           NC = ifelse(is.na(NC), 0, NC),
           G = ifelse(Group == "adult", .5, -.5),
           CT = ifelse(CompType == "New", .5, -.5),
           CM = ifelse(Cond %in% c("CompPresent", "Consolidated",
               "Unconsolidated"), .5, -.5))

dp_item <- dp_frames %>%
    inner_join(tdat2, "TrialID") %>%
    mutate(Comp = ifelse(ID == "Comp", "C", "NC"),
           bin = floor((FID + 1) / 3) * 3) %>%
    count(Item = Sound, Group, CompType, Cond, FID = bin, Comp) %>%
    ungroup() %>%
    spread(Comp, n) %>%
    mutate(C = ifelse(is.na(C), 0, C),
           NC = ifelse(is.na(NC), 0, NC),
           G = ifelse(Group == "adult", .5, -.5),
           CT = (CompType != "New") - mean(CompType == "New"), # bc unbalanced
           CM = ifelse(Cond %in% c("CompPresent", "Consolidated",
               "Unconsolidated"), .5, -.5))

cond_lookup <- dp_subj %>% distinct(Cond) %>%
    mutate(Cond = factor(Cond),
           Condition = c("Control", "Competitor", 
               "Trained on Day 1", "Trained on Day 2", "Untrained"))

group_lookup <- data_frame(Group = factor(c("adult", "child")),
                           Group2 = factor(c("Adults", "Children")))

nperm_runs <- as.integer(commandArgs(TRUE)[1])
nmc <- as.integer(commandArgs(TRUE)[2])

stopifnot(!is.na(nperm_runs) && !is.na(nmc))

item_codes <- dp_item %>%
    select(Item, CompType, CT) %>% distinct()

eff_keep <- c("G", "CT", "CM", "G:CT", "G:CM", "CT:CM", "G:CT:CM")

orig_coef <- daply(dp_item, .(FID), fit_mnom)
orig_pvals <- get_boot_pvals(dp_item, eff_keep, "Item", "CompType", nmc)
saveRDS(orig_coef, file = "derived/results/DP_overall_CR_analysis_orig_by_item.rds")
saveRDS(orig_pvals, file = "derived/results/DP_overall_CR_analysis_orig_pvals_by_item.rds")

ct_eff <- c("CT", "G:CT", "CT:CM", "G:CT:CM")

source("scripts/cluster_setup.R")  # defines 'cl', cores for parallel computing

clusterCall(cl, function(x) {library("plyr"); library("dplyr")}) %>%
  invisible()
clusterExport(cl, setdiff(ls(), "cl"))
cts_list <- parLapply(cl, seq_len(nperm_runs), function(ix) {
              permute_unit(dp_item, item_codes, "CT") %>%
                  get_boot_pvals(ct_eff, "Item", "CompType", nmc) %>%
                  aaply(2, get_clusters, max_cms_only = TRUE)
          })
cts_px <- do.call("rbind", cts_list)

stopCluster(cl)

saveRDS(cts_px, file = "derived/results/DP_overall_CR_analysis_comptype_by_item.rds")
