## clinladder: Clinical Evidence Ladder Visualization
##
## Visual grammar: "Evidence ladders" where each diagnostic or therapeutic
## conclusion is a rung on a ladder, and the rung's visual properties encode:
##   - Evidence strength (strong/moderate/weak)
##   - Evidence type (molecular/imaging/clinical/computational)
##   - Concordance (agree/conflict/missing across modalities)
##   - Actionability (immediately actionable → research-only)
##   - Confidence interval width
##
## Figures generated:
##   A. Oncology evidence ladder: multi-modal evidence for treatment decisions
##      (genomics + pathology + imaging + clinical)
##   B. Rare disease diagnostic odyssey: cascading evidence across tests
##   C. Biomarker validation ladder: evidence accumulation for a panel of biomarkers

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Palettes ─────────────────────────────────────────────────────────────────
clinladder_strength_pal <- function() c(Strong="#264653",Moderate="#2a9d8f",Weak="#e9c46a",Missing="#dee2e6")
clinladder_type_pal <- function() c(Molecular="#264653",Imaging="#2a9d8f",Clinical="#e76f51",
                                     Computational="#7b2cbf",Pathology="#bc6c25")
clinladder_concordance_pal <- function() c(Agree="#2a9d8f",Conflict="#e76f51",Partial="#e9c46a",Missing="#dee2e6")
clinladder_action_pal <- function() c("Tier 1"="#264653","Tier 2"="#2a9d8f","Tier 3"="#e9c46a","Research"="#f4a261")

# ── Data A: Oncology evidence ladder ─────────────────────────────────────────
clinladder_onco_data <- function(seed = 42) {
  set.seed(seed)
  decisions <- c(
    "EGFR TKI eligible","PD-L1 checkpoint eligible","TMB-high classification",
    "HRD-positive status","MSI-high status","KRAS G12C targetable",
    "TP53 LOF (prognosis)","BRCA2 germline (PARP-i)","PIK3CA hotspot (alpelisib)",
    "FGFR2 fusion (erdafitinib)","MET amplification","CDK4/6i eligible",
    "Platinum sensitivity predicted","Immune desert phenotype","Angiogenesis score high")

  modalities <- c("Molecular","Imaging","Clinical","Pathology","Computational")
  nd <- length(decisions); nm <- length(modalities)

  ev_strength <- matrix(sample(c("Strong","Moderate","Weak","Missing"),
                                nd*nm, replace=TRUE, prob=c(.3,.3,.2,.2)),
                         nd, nm, dimnames=list(decisions, modalities))
  # make some decisions have specific patterns
  ev_strength["EGFR TKI eligible",] <- c("Strong","Moderate","Strong","Strong","Moderate")
  ev_strength["PD-L1 checkpoint eligible",] <- c("Strong","Missing","Moderate","Strong","Weak")
  ev_strength["MSI-high status",] <- c("Strong","Missing","Missing","Strong","Moderate")
  ev_strength["Immune desert phenotype",] <- c("Weak","Strong","Moderate","Strong","Strong")

  ev_long <- as.data.frame(as.table(ev_strength), stringsAsFactors=FALSE)
  names(ev_long) <- c("decision","modality","strength")

  # concordance per decision
  conc <- sapply(decisions, function(dec) {
    ss <- ev_strength[dec, ev_strength[dec,] != "Missing"]
    if (length(ss) <= 1) return("Missing")
    if (all(ss %in% c("Strong","Moderate"))) return("Agree")
    if (any(ss == "Strong") && any(ss == "Weak")) return("Conflict")
    return("Partial")
  })

  # actionability tier
  tier <- c("Tier 1","Tier 1","Tier 2","Tier 2","Tier 1","Tier 1","Tier 3","Tier 1",
            "Tier 2","Tier 2","Tier 3","Tier 2","Research","Research","Research")

  # confidence (higher = more strong evidence)
  conf <- sapply(decisions, function(dec) {
    ss <- ev_strength[dec,]
    (.5*sum(ss=="Strong") + .3*sum(ss=="Moderate") + .1*sum(ss=="Weak")) / nm
  })

  list(evidence=ev_long,
       summary=data.frame(decision=decisions, concordance=unname(conc),
                          tier=tier, confidence=unname(conf), stringsAsFactors=FALSE))
}

# ── Data B: Rare disease diagnostic odyssey ──────────────────────────────────
clinladder_rare_data <- function(seed = 7) {
  set.seed(seed)
  tests <- c("Clinical phenotype","Metabolic panel","Karyotype","CMA/SNP array",
             "Gene panel (100 genes)","Exome sequencing","RNA-seq","Trio WGS",
             "Functional assay","Expert consensus")
  n <- length(tests)
  # cumulative evidence score (ramps up then may plateau)
  cum_evidence <- cumsum(c(.15,.08,.05,.12,.18,.25,.10,.30,.15,.20))
  cum_evidence <- cum_evidence / max(cum_evidence)

  # each test has a result
  results <- c("Dysmorphic features","Normal","Normal","16p11.2 del (uncertain)",
               "VUS in SHANK3","Likely pathogenic SHANK3","Aberrant splicing confirmed",
               "De novo confirmed","Protein truncation validated","Phelan-McDermid confirmed")
  significance <- c("Weak","Missing","Missing","Moderate","Moderate","Strong","Strong",
                     "Strong","Strong","Strong")
  time_months <- c(0, 1, 2, 4, 8, 14, 18, 24, 28, 30)

  data.frame(test=factor(tests, levels=tests), result=results,
             significance=significance, time_months=time_months,
             cum_evidence=cum_evidence, stringsAsFactors=FALSE)
}

# ── Data C: Biomarker validation ladder ──────────────────────────────────────
clinladder_biomarker_data <- function(seed = 99) {
  set.seed(seed)
  markers <- c("CD274 (PD-L1)","TMB","MSI","CD8 TIL density","IFN-gamma sig",
               "TGF-beta sig","Angiogenesis","Hypoxia score","CYT score",
               "GEP (18-gene)")
  stages <- c("Discovery","Validation","Clinical trial","Approved/Routine")
  nm <- length(markers); ns <- length(stages)

  # each marker has evidence at each stage
  mat <- matrix(0, nm, ns, dimnames=list(markers, stages))
  mat[,1] <- stats::runif(nm, .3, 1)   # discovery: most have some evidence
  mat[,2] <- mat[,1] * stats::runif(nm, .4, .9)  # validation: some attrition
  mat[,3] <- mat[,2] * stats::runif(nm, .2, .8)  # trial: more attrition
  mat[,4] <- mat[,3] * ifelse(stats::runif(nm) > .5, stats::runif(nm,.3,.7), 0)  # approval: few pass

  # known markers have strong evidence throughout
  mat["CD274 (PD-L1)",] <- c(.95, .88, .82, .78)
  mat["TMB",] <- c(.90, .70, .65, .55)
  mat["MSI",] <- c(.92, .85, .80, .75)

  ev_long <- as.data.frame(as.table(mat), stringsAsFactors=FALSE)
  names(ev_long) <- c("marker","stage","evidence")
  ev_long$stage <- factor(ev_long$stage, levels=stages)
  ev_long$marker <- factor(ev_long$marker, levels=rev(markers))
  ev_long
}

# ── Plot A: Oncology evidence ladder ─────────────────────────────────────────
clinladder_plot_onco <- function() {
  dat <- clinladder_onco_data()
  ev <- dat$evidence; sm <- dat$summary
  str_pal <- clinladder_strength_pal(); conc_pal <- clinladder_concordance_pal()
  act_pal <- clinladder_action_pal()

  # order decisions by confidence
  sm <- sm[order(sm$confidence, decreasing=TRUE),]
  ev$decision <- factor(ev$decision, levels=rev(sm$decision))
  sm$decision <- factor(sm$decision, levels=rev(sm$decision))

  # heatmap: modality × decision, colored by strength
  p_heat <- ggplot(ev, aes_string(x="modality",y="decision",fill="strength")) +
    geom_tile(color="white",linewidth=.5) +
    scale_fill_manual(values=str_pal,name="Evidence\nstrength") +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(angle=35,hjust=1)) +
    labs(x=NULL,y=NULL,title="A1. Multi-modal evidence heatmap")

  # concordance ribbon
  p_conc <- ggplot(sm, aes_string(x="concordance",y="decision",fill="concordance")) +
    geom_tile(color="white",linewidth=.5) +
    scale_fill_manual(values=conc_pal,name="Concordance") +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(),axis.text.y=element_blank()) +
    labs(x=NULL,y=NULL,title="Concordance")

  # actionability tier
  p_tier <- ggplot(sm, aes_string(x="tier",y="decision",fill="tier")) +
    geom_tile(color="white",linewidth=.5) +
    scale_fill_manual(values=act_pal,name="Actionability") +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(),axis.text.y=element_blank()) +
    labs(x=NULL,y=NULL,title="Tier")

  # confidence lollipop
  p_conf <- ggplot(sm, aes_string(x="confidence",y="decision")) +
    geom_segment(aes_string(xend="0",yend="decision"),color="#adb5bd",linewidth=.6) +
    geom_point(aes_string(fill="concordance"),shape=21,size=3,stroke=.4,color="#2b2d42") +
    scale_fill_manual(values=conc_pal,guide="none") +
    geom_vline(xintercept=.5,linetype="dashed",color="#d00000",linewidth=.3) +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank(),axis.text.y=element_blank()) +
    labs(x="Composite confidence",y=NULL,title="Confidence")

  patchwork::wrap_plots(p_heat, p_conc, p_tier, p_conf, nrow=1,
                        widths=c(3, .8, .8, 1.5)) +
    patchwork::plot_annotation(title="clinladder: Oncology evidence ladder for precision treatment",
                               subtitle="Each decision rung shows evidence across modalities, concordance, actionability, and composite confidence")
}

# ── Plot B: Rare disease diagnostic odyssey ──────────────────────────────────
clinladder_plot_rare <- function() {
  d <- clinladder_rare_data()
  sig_pal <- clinladder_strength_pal()

  # timeline: evidence accumulation
  p_timeline <- ggplot(d, aes_string(x="time_months",y="cum_evidence")) +
    geom_area(fill="#264653",alpha=.15) +
    geom_line(color="#264653",linewidth=.8) +
    geom_point(aes_string(fill="significance"),shape=21,size=4,stroke=.5,color="#2b2d42") +
    geom_text(aes_string(label="test"),angle=25,hjust=-.08,vjust=-.5,size=2.5) +
    scale_fill_manual(values=sig_pal,name="Result significance") +
    scale_x_continuous(breaks=d$time_months) +
    scale_y_continuous(limits=c(0,1.15)) +
    theme_minimal(base_size=10) + theme(panel.grid.minor=element_blank()) +
    labs(x="Months from first presentation",y="Cumulative evidence score",
         title="B1. Diagnostic odyssey timeline",
         subtitle="Each point = a test; fill = significance of result; area = cumulative evidence")

  # evidence ladder view: ordered rungs
  d$test_ordered <- factor(d$test, levels=rev(levels(d$test)))
  p_ladder <- ggplot(d, aes_string(y="test_ordered")) +
    geom_segment(aes_string(x="0",xend="cum_evidence",yend="test_ordered",color="significance"),
                 linewidth=3) +
    geom_point(aes_string(x="cum_evidence"),shape=21,size=3,fill="white",stroke=.8,color="#264653") +
    geom_text(aes_string(x="cum_evidence",label="result"),hjust=-.05,size=2.5) +
    scale_color_manual(values=sig_pal,name="Significance") +
    scale_x_continuous(limits=c(0,1.6)) +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank()) +
    labs(x="Cumulative evidence",y=NULL,title="B2. Diagnostic evidence ladder",
         subtitle="Each rung = a test; bar length = evidence contributed; color = significance")

  patchwork::wrap_plots(p_timeline, p_ladder, nrow=2, heights=c(1.5,2)) +
    patchwork::plot_annotation(title="clinladder: Rare disease diagnostic odyssey",
                               subtitle="From phenotype to molecular diagnosis: 30 months, 10 tests, escalating evidence")
}

# ── Plot C: Biomarker validation ladder ──────────────────────────────────────
clinladder_plot_biomarker <- function() {
  d <- clinladder_biomarker_data()

  # heatmap: markers × stages
  p_heat <- ggplot(d, aes_string(x="stage",y="marker",fill="evidence")) +
    geom_tile(color="white",linewidth=.5) +
    scale_fill_gradient2(low="#f8f9fa",mid="#e9c46a",high="#264653",midpoint=.5,
                         limits=c(0,1),name="Evidence\nscore") +
    geom_text(aes_string(label="sprintf('%.2f',evidence)"),size=2.5,color="#2b2d42") +
    theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
                                         axis.text.x=element_text(angle=30,hjust=1)) +
    labs(x=NULL,y=NULL,title="C1. Biomarker validation ladder heatmap")

  # attrition funnel per marker (how much evidence survives each stage)
  surv <- d %>% group_by(marker) %>%
    arrange(stage) %>%
    mutate(pct_retained = evidence / first(evidence))
  surv_wide <- surv %>% filter(stage %in% c("Discovery","Approved/Routine")) %>%
    dplyr::select(marker, stage, evidence) %>%
    pivot_wider(names_from=stage,values_from=evidence)
  surv_wide$attrition <- 1 - surv_wide$`Approved/Routine` / surv_wide$Discovery
  surv_wide <- surv_wide[order(surv_wide$attrition),]
  surv_wide$marker <- factor(surv_wide$marker, levels=surv_wide$marker)

  p_attr <- ggplot(surv_wide, aes_string(x="attrition",y="marker")) +
    geom_segment(aes_string(xend="0",yend="marker"),color="#e76f51",linewidth=2) +
    geom_point(size=3,color="#264653") +
    geom_vline(xintercept=.5,linetype="dashed",color="#d00000",linewidth=.3) +
    scale_x_continuous(labels=scales::percent,limits=c(0,1.05)) +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank()) +
    labs(x="Attrition rate (Discovery → Approved)",y=NULL,
         title="C2. Evidence attrition funnel",
         subtitle="Markers right of dashed line lost >50% of evidence")

  patchwork::wrap_plots(p_heat, p_attr, nrow=1, widths=c(2,1.5)) +
    patchwork::plot_annotation(title="clinladder: Biomarker evidence accumulation and attrition",
                               subtitle="From discovery through clinical validation: which biomarkers survive the evidence gauntlet?")
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
clinladder_plot_demo <- function() { clinladder_plot_onco() }

clinladder_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"clinladder_01_onco.png"),
                  plot=clinladder_plot_onco(), width=15,height=10,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"clinladder_02_rare_disease.png"),
                  plot=clinladder_plot_rare(), width=14,height=11,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"clinladder_03_biomarker.png"),
                  plot=clinladder_plot_biomarker(), width=14,height=9,dpi=300,bg="white")
}
