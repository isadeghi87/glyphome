## counterpath: Temporal Counterfactual Bioinformatics Visualization
##
## Visual grammar: "Decision fans" that show what DID happen alongside what
## COULD HAVE happened if a biological or treatment branch had gone differently.
##
## Three temporal modes:
##   - Observed: solid marks with full confidence
##   - Inferred (model-predicted past): dashed marks with decaying confidence
##   - Counterfactual (alternate timeline): translucent fans that diverge from
##     branch points
##
## Figures generated:
##   A. Treatment response timeline: observed + counterfactual outcome fans for
##      responders vs non-responders
##   B. Gene regulatory cascade: observed cascade + counterfactual "what if
##      transcription factor was silenced" branches
##   C. Cell fate divergence: pseudotime paths with counterfactual fates

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Palettes ─────────────────────────────────────────────────────────────────
counterpath_mode_pal <- function() c(Observed="#264653",Inferred="#2a9d8f",Counterfactual="#e76f51")
counterpath_outcome_pal <- function() c(Responder="#2a9d8f",NonResponder="#e76f51",
                                         "Counter-Responder"="#577590","Counter-NonResp"="#f4a261")

# ── Data A: Treatment response timeline ──────────────────────────────────────
counterpath_treatment_data <- function(seed = 42) {
  set.seed(seed)
  tp <- c(-2, -1, 0, 1, 2, 4, 8, 12, 24)  # weeks
  np <- 6   # patients per group

  mk_patient <- function(pid, group) {
    base <- ifelse(group=="Responder", 85, 80) + stats::rnorm(1, 0, 5)
    slope <- ifelse(group=="Responder", -3.5, .8) + stats::rnorm(1, 0, .5)
    cf_slope <- ifelse(group=="Responder", .5, -2.8) + stats::rnorm(1, 0, .5)
    obs <- pmax(0, base + slope * tp + stats::rnorm(length(tp), 0, 2))
    obs[tp < 0] <- base + stats::rnorm(sum(tp < 0), 0, 2)
    cf <- pmax(0, base + cf_slope * pmax(0, tp) + stats::rnorm(length(tp), 0, 3))
    cf[tp <= 0] <- NA
    conf_obs <- ifelse(tp >= 0, pmax(.3, 1 - .02 * tp), .6 + .15 * tp)
    conf_cf <- ifelse(tp > 0, pmax(.1, .8 - .04 * tp), NA)
    data.frame(patient=pid, group=group, week=tp,
               tumor_burden=obs, counterfactual=cf,
               conf_obs=conf_obs, conf_cf=conf_cf,
               mode=ifelse(tp<0,"Inferred","Observed"), stringsAsFactors=FALSE)
  }
  d <- do.call(rbind, c(
    lapply(1:np, function(i) mk_patient(paste0("R",i), "Responder")),
    lapply(1:np, function(i) mk_patient(paste0("NR",i), "NonResponder"))
  ))
  d
}

# ── Data B: Gene regulatory cascade ─────────────────────────────────────────
counterpath_cascade_data <- function(seed = 7) {
  set.seed(seed)
  steps <- 0:6
  genes <- c("STAT1","IRF1","CXCL10","GBP1","IDO1","PDL1","GZMB")
  # observed cascade: sequential activation
  obs_expr <- matrix(0, length(genes), length(steps), dimnames=list(genes, paste0("t",steps)))
  for (i in seq_along(genes)) {
    onset <- i - 1
    for (j in seq_along(steps)) {
      t <- steps[j]
      if (t >= onset) obs_expr[i,j] <- pmin(1, (t - onset) * .5 + stats::rnorm(1, 0, .05))
      else obs_expr[i,j] <- stats::rnorm(1, .05, .02)
    }
  }
  # counterfactual: STAT1 silenced → downstream attenuated
  cf_expr <- obs_expr
  cf_expr["STAT1",] <- obs_expr["STAT1",] * .1 + stats::rnorm(ncol(cf_expr), 0, .02)
  for (i in 2:length(genes)) cf_expr[i,] <- obs_expr[i,] * (0.15 + .12*(i-2)) + stats::rnorm(ncol(cf_expr), 0, .03)
  cf_expr[] <- pmax(0, pmin(1, cf_expr))

  obs_long <- as.data.frame(as.table(obs_expr), stringsAsFactors=FALSE)
  names(obs_long) <- c("gene","time","expression"); obs_long$mode <- "Observed"
  cf_long <- as.data.frame(as.table(cf_expr), stringsAsFactors=FALSE)
  names(cf_long) <- c("gene","time","expression"); cf_long$mode <- "Counterfactual"
  d <- rbind(obs_long, cf_long)
  d$step <- as.integer(gsub("t","", d$time))
  d$gene <- factor(d$gene, levels=genes)
  d
}

# ── Data C: Cell fate divergence ─────────────────────────────────────────────
counterpath_fate_data <- function(seed = 99) {
  set.seed(seed)
  n <- 100; pt <- seq(0, 1, length.out=n)
  bp <- 45  # branch point index
  # trunk
  tx <- pt[1:bp] * 4; ty <- sin(pt[1:bp] * pi * .6)
  # observed branch A (effector)
  ea_t <- pt[bp:n]; ea_x <- tx[bp] + (ea_t - pt[bp]) * 4; ea_y <- ty[bp] + (ea_t - pt[bp]) * 3
  # observed branch B (exhausted)
  eb_t <- pt[bp:n]; eb_x <- tx[bp] + (eb_t - pt[bp]) * 3.5; eb_y <- ty[bp] - (eb_t - pt[bp]) * 2.5
  # counterfactual: branch A had gone to memory instead
  cf_x <- tx[bp] + (ea_t - pt[bp]) * 5; cf_y <- ty[bp] + (ea_t - pt[bp]) * .5
  confid <- pmax(.1, 1 - 2.5*(ea_t - pt[bp]))

  mk <- function(x,y,t,br,md,cf_v) {
    nn <- length(x)
    data.frame(x=x+stats::rnorm(nn,0,.06), y=y+stats::rnorm(nn,0,.06),
               pseudotime=t, branch=br, mode=md, confidence=cf_v, stringsAsFactors=FALSE)
  }
  rbind(mk(tx,ty,pt[1:bp],"Trunk","Observed",rep(1,bp)),
        mk(ea_x,ea_y,ea_t,"Effector","Observed",rep(1,length(ea_t))),
        mk(eb_x,eb_y,eb_t,"Exhausted","Observed",rep(1,length(eb_t))),
        mk(cf_x,cf_y,ea_t,"Memory (CF)","Counterfactual",confid))
}

# ── Plot A: Treatment response ───────────────────────────────────────────────
counterpath_plot_treatment <- function() {
  d <- counterpath_treatment_data()
  mp <- counterpath_mode_pal()

  # observed + inferred lines per patient
  p_main <- ggplot(d, aes_string(x="week",y="tumor_burden",group="patient")) +
    # counterfactual fans
    geom_ribbon(data=d[!is.na(d$counterfactual),],
                aes_string(ymin="pmin(tumor_burden,counterfactual)",
                           ymax="pmax(tumor_burden,counterfactual)"),
                fill="#e76f51",alpha=.08) +
    geom_line(data=d[!is.na(d$counterfactual),],
              aes_string(y="counterfactual"),linetype="dotted",color="#e76f51",alpha=.4,linewidth=.5) +
    geom_line(data=function(x) x[x$mode=="Observed",], color="#264653",linewidth=.6,alpha=.7) +
    geom_line(data=function(x) x[x$mode=="Inferred",], color="#2a9d8f",linetype="dashed",linewidth=.6,alpha=.7) +
    geom_vline(xintercept=0, linetype="dashed",color="#adb5bd") +
    annotate("text",x=0.5,y=max(d$tumor_burden,na.rm=TRUE)*.95,label="Treatment\nstart",size=2.5,hjust=0) +
    facet_wrap(~group, nrow=1) +
    theme_minimal(base_size=10) +
    theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
          legend.position="top") +
    labs(x="Weeks from treatment",y="Tumor burden (mm²)",
         title="A. Treatment response with counterfactual timelines",
         subtitle="Orange fan: what would have happened under the other group's trajectory")

  # confidence decay ribbon
  conf_d <- d %>% filter(!is.na(conf_cf)) %>% group_by(group,week) %>%
    summarise(mean_cf=mean(conf_cf), se_cf=stats::sd(conf_cf)/sqrt(n()), .groups="drop")
  p_conf <- ggplot(conf_d, aes_string(x="week",y="mean_cf",color="group",fill="group")) +
    geom_ribbon(aes_string(ymin="pmax(0,mean_cf-se_cf)",ymax="pmin(1,mean_cf+se_cf)"),alpha=.15,color=NA) +
    geom_line(linewidth=.8) +
    scale_color_manual(values=c(Responder="#2a9d8f",NonResponder="#e76f51")) +
    scale_fill_manual(values=c(Responder="#2a9d8f",NonResponder="#e76f51")) +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank(),legend.position="none") +
    labs(x="Weeks",y="Counterfactual confidence",title="Confidence decay of counterfactual predictions")

  patchwork::wrap_plots(p_main, p_conf, nrow=2, heights=c(2.5,1)) +
    patchwork::plot_annotation(title="counterpath: Treatment response with temporal counterfactuals",
                               subtitle="Solid = observed; dashed = inferred past; dotted+fan = counterfactual alternate timeline")
}

# ── Plot B: Gene regulatory cascade ─────────────────────────────────────────
counterpath_plot_cascade <- function() {
  d <- counterpath_cascade_data()
  mp <- counterpath_mode_pal()

  # heatmap: observed vs counterfactual side by side
  p_heat <- ggplot(d, aes_string(x="step",y="gene",fill="expression")) +
    geom_tile(color="white",linewidth=.3) +
    facet_wrap(~mode, nrow=1) +
    scale_fill_gradient2(low="#264653",mid="#f8f9fa",high="#e76f51",midpoint=.5,
                         limits=c(0,1),name="Expr") +
    theme_minimal(base_size=10) +
    theme(panel.grid=element_blank(),strip.text=element_text(face="bold")) +
    labs(x="Time step",y=NULL,title="B. Gene cascade: observed vs STAT1-silenced counterfactual")

  # line traces for key genes
  key_genes <- c("STAT1","IRF1","CXCL10","PDL1","GZMB")
  dk <- d[d$gene %in% key_genes,]
  dk$grp <- paste0(dk$gene, "_", dk$mode)
  p_line <- ggplot(dk, aes_string(x="step",y="expression",color="gene",linetype="mode",group="grp")) +
    geom_line(linewidth=.8, alpha=.8) +
    geom_point(size=1.5, alpha=.6) +
    scale_linetype_manual(values=c(Observed="solid",Counterfactual="dashed")) +
    theme_minimal(base_size=9) +
    theme(legend.position="right",panel.grid.minor=element_blank()) +
    labs(x="Time step",y="Expression",color="Gene",linetype="Mode",
         title="Key gene traces: observed (solid) vs counterfactual (dashed)")

  # delta heatmap
  obs_w <- d[d$mode=="Observed", c("gene","step","expression")]
  cf_w <- d[d$mode=="Counterfactual", c("gene","step","expression")]
  delta <- merge(obs_w, cf_w, by=c("gene","step"), suffixes=c("_obs","_cf"))
  delta$diff <- delta$expression_obs - delta$expression_cf
  p_delta <- ggplot(delta, aes_string(x="step",y="gene",fill="diff")) +
    geom_tile(color="white",linewidth=.3) +
    scale_fill_gradient2(low="#264653",mid="#f8f9fa",high="#e76f51",midpoint=0,name="Delta") +
    theme_minimal(base_size=9) + theme(panel.grid=element_blank()) +
    labs(x="Time step",y=NULL,title="Impact of STAT1 silencing (obs − counterfactual)")

  patchwork::wrap_plots(p_heat, patchwork::wrap_plots(p_line, p_delta, nrow=1, widths=c(1.5,1)),
                        nrow=2, heights=c(1.5,1.2)) +
    patchwork::plot_annotation(title="counterpath: Gene regulatory cascade with counterfactual silencing",
                               subtitle="What if the upstream driver had been knocked out? Downstream collapse becomes visible")
}

# ── Plot C: Cell fate divergence ─────────────────────────────────────────────
counterpath_plot_fate <- function() {
  d <- counterpath_fate_data()

  obs_d <- d[d$mode=="Observed",]; cf_d <- d[d$mode=="Counterfactual",]
  pal <- c(Trunk="#264653",Effector="#2a9d8f",Exhausted="#e76f51","Memory (CF)"="#577590")

  p_main <- ggplot() +
    # counterfactual fan
    geom_point(data=cf_d, aes_string(x="x",y="y",alpha="confidence"), color="#577590",size=2) +
    geom_path(data=cf_d[order(cf_d$pseudotime),], aes_string(x="x",y="y"),
              color="#577590",linetype="dashed",linewidth=.8,alpha=.5) +
    # observed
    geom_point(data=obs_d, aes_string(x="x",y="y",color="branch"), size=1.5, alpha=.7) +
    scale_color_manual(values=pal,name="Fate") +
    scale_alpha_continuous(range=c(.08,.5),name="CF Confidence") +
    theme_minimal(base_size=10) + theme(panel.grid=element_blank()) +
    labs(x="Trajectory dim 1",y="Trajectory dim 2",
         title="C. Cell fate divergence with counterfactual memory path",
         subtitle="Dashed translucent = cells that could have become memory instead of effector")

  # confidence along pseudotime for CF
  p_conf <- ggplot(cf_d, aes_string(x="pseudotime",y="confidence")) +
    geom_area(fill="#577590",alpha=.3) + geom_line(color="#577590",linewidth=.6) +
    geom_hline(yintercept=.3,linetype="dashed",color="#d00000") +
    theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank()) +
    labs(x="Pseudotime",y="CF Confidence",title="Counterfactual confidence decay",
         subtitle="Below dashed line: speculative prediction")

  patchwork::wrap_plots(p_main, p_conf, nrow=2, heights=c(2.5,1)) +
    patchwork::plot_annotation(title="counterpath: Branching cell fate with counterfactual alternatives",
                               subtitle="What if cells had taken the other path? Confidence fan narrows with temporal distance")
}

# ── Demo + Save ──────────────────────────────────────────────────────────────
counterpath_plot_demo <- function() { counterpath_plot_treatment() }

counterpath_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  ggplot2::ggsave(file.path(output_dir,"counterpath_01_treatment.png"),
                  plot=counterpath_plot_treatment(), width=14,height=10,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"counterpath_02_cascade.png"),
                  plot=counterpath_plot_cascade(), width=14,height=10,dpi=300,bg="white")
  ggplot2::ggsave(file.path(output_dir,"counterpath_03_fate.png"),
                  plot=counterpath_plot_fate(), width=12,height=10,dpi=300,bg="white")
}
