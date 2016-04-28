library("plyr")
library("dplyr")
library("tidyr")
library("Rcpp")

## pull in information about each trial
extractTrialInfo <- function(x) {
    ff <- read.delim(as.character(x$fname), stringsAsFactors = FALSE)
    ff2 <- ff %>% mutate(Resp = as.numeric(sub("^AOI", "", CRESP)),
                         TLoc = as.numeric(sub("^AOI", "", CorrectAnswer)),
                         Acc = ACC) %>%
           select(Subject, Sound, AOI1, AOI2, AOI3, AOI4, Condition,
                  TLoc, Resp, Acc, RT) %>%
           unique()
    ff2$TOrd <- 1:nrow(ff2)
    ff2
}

## make gaze data cumulative to target
makeCumulative <- function(x, maxWin) {
    if (max(x$FID) == maxWin) {
        res <- x 
        res$Pad <- 0
    } else {
        x$Pad <- 0
        res <- rbind(x,
                     data.frame(TrialID = x$TrialID[1],
                                FID = seq(max(x$FID) + 1, maxWin, 1),
                                ID = "Targ",
                                Pad = 1))
    }
    return(res)
}

## pull in the eyedata
extractEyeData <- function(x) {
    ff <- read.delim(as.character(x$fname), stringsAsFactors = FALSE) %>%
        ## only Preview/StimSlide
        filter(UserDefined_1 != "Fixation" & UserDefined_1 != "") %>%
        ## change "SquareX" to  number X        
        mutate(AOI = as.integer(sub("^Square","", AOI)), 
               Val = ValidityRightEye * ValidityLeftEye) %>%
        select(Subject, Sound, FID = ID, AOI, Phase = UserDefined_1,
               x = CursorX, y = CursorY) %>%
        inner_join(onsetTimes, by = c("Sound")) %>%
        ## set zero-frame to onset of the word? (article?)
        mutate(FOnset = round(60 * Onset / 1000)) ## <-- ** CHANGE THIS!

    ff2 <- ff %>%
        group_by(Subject, Sound) %>%
        filter(good = "StimSlide" %in% Phase) %>%
        mutate(t0 = min(ifelse(Phase == "StimSlide",
                   FID + FOnset, NA), na.rm = TRUE))

    if (ff2 %>%
            select(Sound, t0) %>%
            unique() %>%
            group_by(Sound) %>%
            filter(n() != 1) %>%
            nrow(.) != 0) {
        stop("couldn't find soundfile for subject in row ", x$drow)
    } else {}
    # remap positions based on Eye data
    ff2$MapAOI <- inRegion(aoipos, ff2$x, ff2$y)
    # now map to AOI
    ff3 <- ff2 %>%
        mutate(FID = FID - t0) %>% ungroup() %>% 
        inner_join(select(tdat2, TrialID, Subject, Sound),
                   by=c("Subject", "Sound")) %>%
        select(TrialID, FID, AOI, MapAOI, x, y) %>%
        left_join(trialaoi, by = c(TrialID = "TrialID", MapAOI = "Loc")) %>%
        rename(xyID = Role) %>% select(-Image)
    # now map eye to AOI from the Tobii output
    ff4 <- ff3 %>% 
        left_join(trialaoi, by = c(TrialID = "TrialID", AOI = "Loc")) %>%
        rename(tobiiID = Role) %>%
        select(TrialID, FID, xyID, tobiiID) %>%
        arrange(TrialID, FID)
    # xyID: aoi identity mapped from x/y coords
    # tobiiID: aoi identity mapped from tobii output
    ff4$xyID <- ff4$xyID
    ff4$tobiiID <- ff4$tobiiID
    return(ff4)
}

## pull in trial information from the eyedata files
adults <- list.files("gazedata/adults", full.names = TRUE)
kids <- list.files("gazedata/children", full.names = TRUE)
dfiles <- data.frame(Group = rep(c("adult", "child"),
                         c(length(adults), length(kids))),
                     fname = c(adults, kids),
                     drow = 1:(length(adults) + length(kids)))

## load in file with information about each of the displays
stim <- read.csv("FIXED_stimulus_list.csv", stringsAsFactors = FALSE) %>%
    filter(Condition != "novelfiller") %>%
    select(TLoc = CorrectAnswer, Sound = Target, Condition,
           Targ = Target.AOI, Comp = Competitor.AOI, DN = NovelDistractor,
           DE = ExistingDistrator, AOI1 = Image1, AOI2 = Image2, AOI3 = Image3,
           AOI4 = Image4, Onset = TargetOnset) %>%
    unique() %>% arrange(Sound, Condition) %>%
    mutate(DispID = row_number(),
           AOI1 = as.character(AOI1),
           AOI2 = as.character(AOI2),
           AOI3 = as.character(AOI3),
           AOI4 = as.character(AOI4))

onsetTimes <- read.csv("onsetTimes.csv", stringsAsFactors = FALSE)
## ## this is the code that was originally used to generate the 
## ## onsetTimesTemplate:
## onsetTimes <- stim %>% filter(Condition!="novelfiller") %>%
##     group_by(Sound) %>% summarize(Onset=round(mean(Onset))) %>%
##    mutate(Sound=as.character(Sound))

aoiloc0 <- ddply(stim, .(DispID, Sound), function(x) {
    toCol <- function(x) {paste0("AOI", x)}
    data.frame(Role = c("Targ", "Comp", "DN", "DE"),
               Loc = c(x$Targ, x$Comp, x$DN, x$DE),
               Image = c(x[, toCol(x$Targ)], x[, toCol(x$Comp)],
                   x[, toCol(x$DN)], x[, toCol(x$DE)]))
})
aoiloc <- aoiloc0 %>% mutate(Role = as.character(Role),
                             Image = as.character(Image))


## check for errors
if (aoiloc %>% select(Sound, Image, Role) %>% unique() %>%
    group_by(Sound, Image) %>%
    filter(n() > 1) %>% arrange(Sound, Image) %>% nrow(.) == 0) {
} else {
    stop("error: problem with 'aoiloc'?")
}

aoirole <- aoiloc %>% select(Sound, Image, Role) %>% unique() %>%
    arrange(Sound, Image, Role)

## extract the trial information from the gaze files
## ignoring the eye/mouse data for now
tdat0 <- dfiles %>%
    group_by(drow, Group) %>%
    do(extractTrialInfo(.)) %>%
    ungroup()

tdat0$TrialID <- 1:nrow(tdat0)
tdat01 <- tdat0 %>%
    group_by(Subject) %>%
    mutate(TOrd = row_number()) %>%
    ungroup()

trialaoi <- tdat01 %>%
    filter(Condition != "novelfiller", Condition !="practice") %>%
    select(TrialID, Sound, AOI1:AOI4) %>%
    gather(Loc, Image, AOI1, AOI2, AOI3, AOI4) %>%
    mutate(Loc = as.numeric(sub("^AOI", "", Loc))) %>% 
    left_join(aoirole, by = c("Sound","Image")) %>%
    arrange(TrialID, Loc) %>%
    ## mutate(Role = factor(Role, levels = c("Targ","Comp","DN","DE"))) %>%
    select(-Sound)
trialaoi$Loc <- as.integer(trialaoi$Loc)

## check it
## test 1: 4 pictures per trialID?
if (trialaoi %>% group_by(TrialID) %>% filter(n() != 4) %>% nrow(.) != 0) {
    stop("should be 4 pictures per trialID")
} else {}
## test 2: correspondence of TLoc in tdat01 to target location in trialaoi?
if (tdat01 %>% select(TrialID, Loc = TLoc) %>% inner_join(trialaoi) %>%
    filter(Role != "Targ") %>% nrow(.) != 0) {
    stop("TLoc didn't correspond to target location in trialaoi")
} else {}
## test 3: each trial has pictures in all four roles
t3 <- ddply(trialaoi, .(TrialID), function(x) {
    paste(sort(x$Role), collapse = "_")
})
if (length(unique(t3$V1)) != 1) {
    stop("each trial should have four rows in 'trialaoi'")
} else {}
rm(t3)

tdat <- tdat01 %>%
    filter(Condition != "novelfiller", Condition != "practice") %>%
    select(Subject, Group, TrialID, TOrd, Sound, Condition, Loc = Resp, RT) %>%
    mutate(Loc = as.integer(Loc)) %>%
    inner_join(trialaoi, by = c("TrialID","Loc")) %>%
    select(Subject, Group, TrialID, TOrd, Sound, Condition, Resp = Role, RT) %>%
    inner_join(onsetTimes, by="Sound") %>%
    mutate(RTfromOnset = RT - Onset,
           DPlag = DP - Onset) %>%
    select(-Article, -RT, -Onset)

missing <- anti_join(tdat, onsetTimes, by = "Sound") %>%
    select(Sound) %>% unique() %>% arrange(Sound)
if (nrow(missing) != 0) {
    stop("some of the onset times seem to be missing")
} else {}

compType <- data.frame(Condition = sort(unique(tdat$Condition)),
                       CompType = c("Existing","New","Existing","New","New"),
                       Cond = c("CompAbsent", "Consolidated","CompPresent",
                           "Unconsolidated","Untrained"),
                       stringsAsFactors = FALSE)

tdat2 <- inner_join(tdat, compType, by="Condition") %>%
    select(-Condition)

sourceCpp("inregion.cpp")

aoipos <- data.frame(Loc2 = 1:4,
                     x1 = c(0, 704, 0, 704), 
                     y1 = c(0, 0, 564, 564),
                     x2 = c(575, 1279, 575, 1279),
                     y2 = c(460, 460, 1023, 1023))

edat <- dfiles %>%
    group_by(drow, Group) %>%
    do(extractEyeData(.)) %>%
    ungroup() %>%
    select(-drow, -Group)

trial.end <- edat %>% 
    group_by(TrialID) %>%
    summarize(maxFr = max(FID), maxMs = maxFr * 1000 / 60) %>%
    inner_join(select(tdat2, TrialID, RTfromOnset), by = "TrialID") %>%
    mutate(diff = maxMs - RTfromOnset)

## hist(trial.end$diff)
## so RTfromOnset precedes end of eye data recording by about 200--450 ms
## so would it make sense to stop recording at RTfromOnset?

tend2 <- trial.end %>% mutate(RTFR = round(60 * RTfromOnset / 1000))

## check which one has more target gaze and use that
ed.test <- tend2 %>% 
    select(TrialID, maxFr, RTFR) %>%
    inner_join(edat, by = "TrialID")

ed.test %>% 
    filter(FID == RTFR) %>% 
    group_by(xyID) %>% 
    summarize(n = n()) %>% ungroup() %>%
    mutate(p = n / sum(n))
ed.test %>% 
    filter(FID == maxFr) %>% 
    group_by(xyID) %>% 
    summarize(n = n()) %>% ungroup() %>%
    mutate(p = n / sum(n))

## OK, so we will cut it off at maxFr
edat2 <- ed.test %>% filter(FID < maxFr) %>% select(TrialID, FID, ID = xyID)
trials <- edat2 %>% select(TrialID) %>% unique() %>% nrow()

## have a look at median response times
tdat2 %>% group_by(Group, CompType, Cond) %>% 
    summarize(medRT = median(RTfromOnset))

## don't go beyond 2.5 seconds (longest median RT for kids;
## Unconsolidated condition)
edat3 <- edat2 %>% filter(FID>=-30 & FID<=150)

## find trials with extremely fast responses <400 ms
edat3 %>% group_by(TrialID) %>% filter(max(FID) < 24) %>%
    ungroup() %>% select(TrialID, FID) %>%
    group_by(TrialID) %>% summarize(max(FID))

## get rid of them
edat4 <- edat3 %>% group_by(TrialID) %>%
    filter(max(FID) >= 24) %>% ungroup()

## make trials cumulative to selection
edat.ctt <- edat4 %>% group_by(TrialID) %>%
    do(makeCumulative(., maxWin=150)) %>%
    ungroup() %>%
    mutate(ID = ifelse(is.na(ID), "X", ID))

saveRDS(tdat2, file = "derived/trial_data.rds")
saveRDS(edat4, file = "derived/eye_data.rds")
saveRDS(edat.ctt, file = "derived/eye_data_cumulative.rds")
saveRDS(onsetTimes, file = "derived/onset_times.rds")
saveRDS(trialaoi, file = "derived/trial_aoi.rds")
