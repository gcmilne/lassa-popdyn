#----------------------------------------------------------------------------------#
# Integral projection model of Mastomys natalensis population & infection dynamics #
#----------------------------------------------------------------------------------#

# FUNCTION TO RUN IPM ----
RunRodentModel <- function(model, env, cmr.data, demog.pars, inf.pars, 
                           n_init = 1000, start_date = "1996-01-01", 
                           end_date = NULL, seed, dt = 28,
                           summarised = TRUE) {
  
  with(as.list(c(demog.pars, inf.pars)),{
    
    # 0. Function argument definitions ----
    
    # model: stat model estimates of climate-demographic parameters for survival, recruitment, growth, & male population size
    # env: decomposed temperature & precipitation data in long format by lag combination 
    #      (columns: date, tseas, tvar, pseas, pvar, t_lag, p_lag)
    # cmr.data: rodent capture-mark-capture data (colnames = )
    # demog.pars: list of demographic parameters with names as follows:
    #           - effective_litter: expected no. (female) offspring in a litter
    #           - survival_adj: recapture modifier
    #           - pop_adj: used to relate model population size to data population size
    #           - juvenile_wgt_range: weight range (in grams) of juveniles
    # inf.pars: list of infection parameters with names as follows:
    #           - inf_prev: initial infection prevalence
    #           - vertical_trans: rate of vertical transmission 
    #           - trans_rate: rate of horizontal transmission
    #           - dens_cont_toggle: turn density-dependent horizontal transmission on (1) of off (0)
    #           - mat_ab_toggle: turn passing of maternal antibodies from recovered mothers to newborns on (1) of off (0)
    # n_init: initial population size
    # start_date: simulation start date in YYYY-MM-DD format
    # end_date simulation end date in YYYY-MM-DD format
    # seed: numeric seed value to ensure repeatability using `sample()`
    # dt: time step (in days) - ** keep fixed at default of 28 days **
    
    # 1. Set up model & initialise variables ----
     
    # Set seed to ensure later use of `sample()` is repeatable
    set.seed(seed)
    
    # Generate a vector of dates to simulate over
    dates <- unique(env$date)
    
    # Ensure all climate data dates are in correct dt
    if (all(diff(dates) != dt)) stop("1≥ climate data date interval not consistent with model dt")
    
    # Define number of climate lag combinations in input climate data
    L <- nrow(env) / length(dates)
    
    # Time step indices
    TT <- 1:length(dates)
    
    # Generate body weight mesh to simulate model over
    dw <- 1  # weight step size (in grams)
    WW <- seq(0, 100, by = dw)

    # Population size adjustment parameter
    pop_adj <- pop_adj * mean(cmr.data$POP) / (n_init * 2)  # n_init*2 to account for males in denominator
    
    # Get body size distribution of offspring: expected no. offspring of each weight
    offspring_dist <- effective_litter * (WW >= juvenile_wgt_range[1] & WW <= juvenile_wgt_range[2]) / 
      length(which(WW >= juvenile_wgt_range[1] & WW <= juvenile_wgt_range[2]))
    
    # 2. Find all the statistical model estimated parameters ----
    
    # Replace NA with 0 if needed
    na_to_zero <- function(x) {
      val <- unique(x)
      if (is.na(val)) 0 else val
    }
    
    # Simple unique grab
    get_unique <- function(x) unique(x)
    
    # Direct assignment (for things like lag_wgt_val)
    get_raw <- function(x) x
    
    # Define variable mappings for each process
    param_map <- list(
      surv = list(
        int    = list(col = "alpha",       fn = get_unique),
        w      = list(col = "beta_weight", fn = get_unique),
        pop    = list(col = "beta_pop",    fn = get_unique),
        ttrend = list(col = "beta_ttrend",  fn = na_to_zero),
        tseas  = list(col = "beta_tseas",  fn = na_to_zero),
        tvar   = list(col = "beta_tvar",   fn = na_to_zero),
        pseas  = list(col = "beta_pseas",  fn = na_to_zero),
        pvar   = list(col = "beta_pvar",   fn = na_to_zero),
        ptrend = list(col = "beta_ptrend",  fn = na_to_zero),
        lagwgt = list(col = "lag_wgt_val", fn = get_raw)
      ),
      rec = list(
        int    = list(col = "alpha",       fn = get_unique),
        w      = list(col = "beta_weight", fn = get_unique),
        pop    = list(col = "beta_pop",    fn = get_unique),
        ttrend = list(col = "beta_ttrend",  fn = na_to_zero),
        tseas  = list(col = "beta_tseas",  fn = na_to_zero),
        tvar   = list(col = "beta_tvar",   fn = na_to_zero),
        pseas  = list(col = "beta_pseas",  fn = na_to_zero),
        pvar   = list(col = "beta_pvar",   fn = na_to_zero),
        ptrend = list(col = "beta_ptrend",  fn = na_to_zero),
        lagwgt = list(col = "lag_wgt_val", fn = get_raw)
      ),
      grow = list(
        int.mu    = list(col = "alpha_mu",    fn = get_unique),
        w.mu      = list(col = "beta_weight", fn = get_unique),
        pop.mu    = list(col = "beta_pop",    fn = get_unique),
        ttrend.mu = list(col = "beta_ttrend",  fn = na_to_zero),
        tseas.mu  = list(col = "beta_tseas",  fn = na_to_zero),
        tvar.mu   = list(col = "beta_tvar",   fn = na_to_zero),
        pseas.mu  = list(col = "beta_pseas",  fn = na_to_zero),
        pvar.mu   = list(col = "beta_pvar",   fn = na_to_zero),
        ptrend.mu = list(col = "beta_ptrend",  fn = na_to_zero),
        int.sigma = list(col = "alpha_sigma", fn = get_unique),
        w.sigma   = list(col = "beta_sigma",  fn = get_unique),
        lagwgt    = list(col = "lag_wgt_val", fn = get_raw)
      )
    )
    
    # Loop over each demographic process and store objects in environment
    for (proc in names(param_map)) {
      proc_pars <- model %>% filter(process == proc)
      for (var in names(param_map[[proc]])) {
        colname <- param_map[[proc]][[var]]$col
        fn <- param_map[[proc]][[var]]$fn
        val <- fn(proc_pars[[colname]])
        assign(paste0(var, ".", proc), val, envir = .GlobalEnv)
      }
    }
    
    ## Male population size
    mp.pars <- model %>% filter(process == "male_pop")
    
    # 3. Get realistic distribution of initial weights ----
    
    # Omit pregnant & lactating rodents to not bias weights
    NP <- cmr.data[!(cmr.data$PRE | cmr.data$LAC), ]
    
    # Sample row numbers
    init <- sample(1:nrow(NP), size = (1 - inf_prev) * n_init, replace = TRUE)
    
    # Use these row numbers to get sample of weights
    init <- NP[init, "WEIGHT"]
    
    if (dw %% 1 != 0) {
      # Add random noise to the weights
      init <- init + sample(seq(-0.5, 0.5, by = dw)[-1], size = n_init, replace = TRUE)
    }
    
    # Generate random sample of infected individuals of various weights
    inf_init <- sample(init, size = inf_prev * n_init, replace = FALSE)
    
    # Calculate initial no. individuals in each weight class
    dens_s_init <- dens(vector = init, values = WW) # adult S
    dens_i_init <- dens(vector = inf_init, values = WW) # adult I
    
    dens_s_init_sub <- dens_s_init[1:(juvenile_wgt_range[2] + 1)] # subadult S
    dens_s_init_sub <- c(dens_s_init_sub, rep(0, length(WW) - length(dens_s_init_sub)))
    dens_i_init_sub <- dens_i_init[1:(juvenile_wgt_range[2] + 1)] # subadult I
    dens_i_init_sub <- c(dens_i_init_sub, rep(0, length(WW) - length(dens_i_init_sub)))
    
    dens_s_init[1:(juvenile_wgt_range[2] + 1)] <- 0
    dens_i_init[1:(juvenile_wgt_range[2] + 1)] <- 0
    
    # 4. Initialise model storage & store model output for time t = 1 ----
    D <- data.frame(time   = 1,
                    date   = start_date,
                    weight = WW,
                    # Track infected, uninfected, & maternal Ab+ juveniles
                    juveniles_i = 0,
                    juveniles_s = 0,
                    # Distribution of weights across all compartments
                    D_all  = dens(vector = init, values = WW),
                    # Pregnant adults (state J)
                    J_s    = 0,
                    J_h    = 0,
                    J_v    = 0,
                    J_r    = 0,
                    J_all  = 0,
                    # Non-pregnant adults
                    D_adult_s   = dens_s_init,
                    D_adult_h   = dens_i_init,
                    D_adult_v   = 0,
                    D_adult_r   = 0,
                    D_adult_all = dens_s_init + dens_i_init,
                    # Sub-adults (2 time steps)
                    Dsub_m   = 0,
                    Dsub_s1  = (1/2) * dens_s_init_sub,
                    Dsub_s2  = (1/2) * dens_s_init_sub,
                    Dsub_h1  = (1/2) * dens_i_init_sub,
                    Dsub_h2  = (1/2) * dens_i_init_sub,
                    Dsub_v1  = 0,
                    Dsub_v2  = 0,
                    Dsub_r1  = 0,
                    Dsub_r2  = 0,
                    Dsub_all = 0,
                    # Survival & pregnancy probabilities, force of infection
                    surv_prob = NA,
                    preg_prob = NA,
                    foi = NA)
    
    # 5. Time loop to calculate changes in demographic & infection states ----
    for (t in TT[-1]) {
      
      date <- (start_date - dt) + dt * t  # Current simulation date
      prev_t <- t - 1                     # Previous time step
      
      ## i. Retrieve population distributions from previous time step ----
      D_prev <- D[D$time == prev_t, ]
      
      # Extract previous states
      D_all_prev <- D_prev$D_all
      
      D_adult_s_prev <- D_prev$D_adult_s
      D_adult_h_prev <- D_prev$D_adult_h
      D_adult_v_prev <- D_prev$D_adult_v
      D_adult_r_prev <- D_prev$D_adult_r

      J_s_prev <- D_prev$J_s
      J_h_prev <- D_prev$J_h
      J_v_prev <- D_prev$J_v
      J_r_prev <- D_prev$J_r
      J_all_prev <- J_s_prev + J_h_prev + J_v_prev + J_r_prev
      
      Dsub_m_prev  <- D_prev$Dsub_m
      Dsub_s1_prev <- D_prev$Dsub_s1
      Dsub_s2_prev <- D_prev$Dsub_s2
      Dsub_h1_prev <- D_prev$Dsub_h1
      Dsub_h2_prev <- D_prev$Dsub_h2
      Dsub_v1_prev <- D_prev$Dsub_v1
      Dsub_v2_prev <- D_prev$Dsub_v2
      Dsub_r1_prev <- D_prev$Dsub_r1
      Dsub_r2_prev <- D_prev$Dsub_r2
      
      Dsub_s_prev  <- Dsub_s1_prev + Dsub_s2_prev
      Dsub_h_prev  <- Dsub_h1_prev + Dsub_h2_prev
      Dsub_v_prev  <- Dsub_v1_prev + Dsub_v2_prev
      Dsub_r_prev  <- Dsub_r1_prev + Dsub_r2_prev

      ## i. Climate data lookup ----
      d.ind <- c(t, t + length(dates) * 1:(L-1))
      ttrend.vec <- env$ttrend[d.ind]
      tseas.vec  <- env$tseas[d.ind]
      tvar.vec   <- env$tvar[d.ind]
      ptrend.vec <- env$ptrend[d.ind]
      pseas.vec  <- env$pseas[d.ind]
      pvar.vec   <- env$pvar[d.ind]
      
      ## ii. Demographic process model calculations ----
      
      # Male population size and total population
      set.seed(seed)
      male_pop <- rtruncnorm(
        n = 1, a = 0, b = Inf, 
        mean = as.numeric(mp.pars$alpha + mp.pars[, paste0("alpha_month_", month(date))] + mp.pars$beta * sum(D_all_prev)),
        sd = mp.pars$sigma
      )
      tot_pop <- sum(D_all_prev) + male_pop
      
      # Survival calculation
      surv.lodds <- int.surv + 
        w.surv      * WW +
        pop.surv    * tot_pop * pop_adj +
        ttrend.surv * sum(ttrend.vec * lagwgt.surv) +
        tseas.surv  * sum(tseas.vec  * lagwgt.surv) +
        tvar.surv   * sum(tvar.vec   * lagwgt.surv) +
        ptrend.surv * sum(ptrend.vec * lagwgt.surv) +
        pseas.surv  * sum(pseas.vec  * lagwgt.surv) +
        pvar.surv   * sum(pvar.vec   * lagwgt.surv)
      
      surv.prob <- 1 / (1 + exp(-surv.lodds))
      S_calc <- survival_adj + (1 - survival_adj) * surv.prob

      # Recruitment calculation
      rec.lodds <- int.rec +
        w.rec      * WW +
        pop.rec    * tot_pop * pop_adj +
        ttrend.rec * sum(ttrend.vec * lagwgt.rec) +
        tseas.rec  * sum(tseas.vec  * lagwgt.rec) +
        tvar.rec   * sum(tvar.vec   * lagwgt.rec) +
        ptrend.rec * sum(ptrend.vec * lagwgt.rec) +
        pseas.rec  * sum(pseas.vec  * lagwgt.rec) +
        pvar.rec   * sum(pvar.vec   * lagwgt.rec)
      
      R_calc <- 1 / (1 + exp(-rec.lodds))

      # Growth distribution
      G_mean <- int.mu.grow + 
        w.mu.grow      * WW +
        pop.mu.grow    * tot_pop * pop_adj +
        ttrend.mu.grow * sum(ttrend.vec * lagwgt.grow) +
        tseas.mu.grow  * sum(tseas.vec  * lagwgt.grow) +
        tvar.mu.grow   * sum(tvar.vec   * lagwgt.grow) +
        ptrend.mu.grow * sum(ptrend.vec * lagwgt.grow) +
        pseas.mu.grow  * sum(pseas.vec  * lagwgt.grow) +
        pvar.mu.grow   * sum(pvar.vec   * lagwgt.grow)
      
      G_sd <- exp(int.sigma.grow + w.sigma.grow * WW)
      G_dist <- dnorm(t(array(WW, c(length(WW), length(WW)))), 
                      mean = array(G_mean, c(length(WW), length(WW))), 
                      sd = array(G_sd, c(length(WW), length(WW))))
      G_dist <- G_dist / rowSums(G_dist)
      
      ## iii. Force of infection ----
      non_inf_size <- sum(D_adult_s_prev + Dsub_s_prev + J_s_prev)
      inf_size     <- sum(D_adult_h_prev + D_adult_v_prev + Dsub_h_prev + Dsub_v_prev + J_h_prev + J_v_prev)
      h_inf <- GetFoI(non_inf_size = non_inf_size, inf_size = inf_size, beta = trans_rate, dens_cont = dens_cont_toggle)
      
      ## iv. Demographic & infection transitions ----
      
      # Pregnant individuals by disease state
      J_s <- as.vector(pmax(0, R_calc * (1 - h_inf) * S_calc * D_adult_s_prev) %*% G_dist)
      J_h <- as.vector(pmax(0, R_calc *      h_inf  * S_calc * D_adult_s_prev) %*% G_dist)
      J_v <- as.vector(pmax(0, R_calc               * S_calc * D_adult_v_prev) %*% G_dist)
      J_r <- as.vector(pmax(0, R_calc               * S_calc * (D_adult_r_prev + D_adult_h_prev)) %*% G_dist)
      
      J_all <- J_s + J_h + J_v + J_r
      
      # Juvenile births
      if (mat_ab_toggle == 1) {  # maternal Ab from recovered pregnancies
        juveniles_s   <- offspring_dist * sum((J_s_prev) * S_calc)
        juveniles_m   <- offspring_dist * sum((J_r_prev) * S_calc)
      } else {  # no maternal antibodies
        juveniles_s   <- offspring_dist * sum((J_s_prev + J_r_prev) * S_calc)
        juveniles_m   <- 0
      }
      juveniles_i   <- offspring_dist * sum((J_h_prev + J_v_prev) * S_calc)
      juveniles_all <- juveniles_s + juveniles_m + juveniles_i
      
      # New sub-adults, & sub-adult demographic & weight transitions
      Dsub_m <- juveniles_m  # Maternal Ab when born from R mother (transient for 1 time step) [DOI ref: 10.1111/1365-2656.13107]
      
      Dsub_s1 <- juveniles_s + (1 - vertical_trans) * juveniles_i # new susceptibles
      Dsub_s2 <- as.vector(pmax(0, (1 - h_inf) * S_calc * Dsub_s1_prev + S_calc * Dsub_m_prev) %*% G_dist)

      Dsub_h1 <- 0  # assume no horizontally infected just after birth
      Dsub_h2 <- as.vector(pmax(0, h_inf * S_calc * Dsub_s1_prev) %*% G_dist)

      Dsub_v1 <- vertical_trans * juveniles_i # new vertically infected
      Dsub_v2 <- as.vector(pmax(0, S_calc * Dsub_v1_prev) %*% G_dist)

      Dsub_r1 <- 0
      Dsub_r2 <- as.vector(pmax(0, S_calc * (Dsub_r1_prev + Dsub_h1_prev)) %*% G_dist)

      Dsub_all <- Dsub_m + Dsub_s1 + Dsub_s2 + Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2 + Dsub_r1 + Dsub_r2
      
      # Update adults & add newly matured sub-adults
      D_adult_s <- as.vector(pmax(0, (1 - h_inf) * S_calc * (D_adult_s_prev + Dsub_s2_prev)) %*% G_dist)
      D_adult_h <- as.vector(pmax(0, h_inf * S_calc * (D_adult_s_prev + Dsub_s2_prev)) %*% G_dist)
      D_adult_v <- as.vector(pmax(0, S_calc * (D_adult_v_prev + Dsub_v2_prev)) %*% G_dist)
      D_adult_r <- as.vector(pmax(0, S_calc * (D_adult_r_prev + D_adult_h_prev + Dsub_r2_prev + Dsub_h2_prev)) %*% G_dist)
      
      D_adult_all <- D_adult_s + D_adult_h + D_adult_v + D_adult_r
      
      # Overall population density
      D_all <- J_all + D_adult_all + Dsub_all
      
      ## v. Store output ----
      D_t <- data.frame(
        time   = t,
        date   = date,
        weight = WW,
        juveniles_i = juveniles_i,
        juveniles_s = juveniles_s,
        D_all  = D_all,
        J_s    = J_s,
        J_h    = J_h,
        J_v    = J_v,
        J_r    = J_r,
        J_all  = J_all,
        D_adult_s = D_adult_s,
        D_adult_h = D_adult_h,
        D_adult_v = D_adult_v,
        D_adult_r = D_adult_r,
        D_adult_all = D_adult_all,
        Dsub_m   = Dsub_m,
        Dsub_s1  = Dsub_s1,
        Dsub_s2  = Dsub_s2,
        Dsub_h1  = Dsub_h1,
        Dsub_h2  = Dsub_h2,
        Dsub_v1  = Dsub_v1,
        Dsub_v2  = Dsub_v2,
        Dsub_r1  = Dsub_r1,
        Dsub_r2  = Dsub_r2,
        Dsub_all = Dsub_all,
        surv_prob = S_calc,
        preg_prob = R_calc, 
        foi = h_inf
      )
      
      # Append to overall output
      D <- rbind(D, D_t)
    }
    
    # Calculate modelled population size, prevalence & seroprevalence
    if (summarised) {
      D <- D %>%
        group_by(time,date) %>%
        summarise(.groups="keep",
          Dsub_m = sum(Dsub_m),
          Dsub_s1 = sum(Dsub_s1),
          Dsub_s2 = sum(Dsub_s2),
          Dsub_h1 = sum(Dsub_h1),
          Dsub_h2 = sum(Dsub_h2),
          Dsub_v1 = sum(Dsub_v1),
          Dsub_v2 = sum(Dsub_v2),
          Dsub_r1 = sum(Dsub_r1),
          Dsub_r2 = sum(Dsub_r2),
          D_adult_s = sum(D_adult_s),
          D_adult_h = sum(D_adult_h),
          D_adult_v = sum(D_adult_v),
          D_adult_r = sum(D_adult_r),
          J_s = sum(J_s),
          J_h = sum(J_h),
          J_v = sum(J_v),
          J_r = sum(J_r),
          J_all = sum(J_all),
          Dsub_all = sum(Dsub_all),
          mean_adult_wgt = sum(D_adult_all*weight)/sum(D_adult_all),
          D_adult_all = sum(D_adult_all),
          D_all = sum(D_all),
          surv_prob = mean(surv_prob),
          preg_prob = mean(preg_prob),
          foi = mean(foi)
        ) %>%
        mutate(
          mod.pop   = D_all,
          sub.pop   = Dsub_all,
          adult.pop = D_adult_all + J_all,
          preg.pop  = J_all,
          nonpreg.pop = D_adult_all,
          inf.sub = (Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2),
          inf.adult = (D_adult_h + D_adult_v + J_h + J_v),
          prev.sub = inf.sub / sub.pop,
          prev.adult = inf.adult / adult.pop,
          pos.sub = (Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2 + Dsub_r1 + Dsub_r2 + Dsub_m),
          pos.adult = (D_adult_h + D_adult_v + D_adult_r + J_h + J_v + J_r),
          seroprev.sub = pos.sub / sub.pop,
          seroprev.adult = pos.adult / adult.pop,
          births.total = Dsub_m + Dsub_s1 + Dsub_h1 + Dsub_v1 + Dsub_r1,
          deaths.total = 0,
          infectious.laststep = 0
        ) %>%
        ungroup()
      
      # deaths since the last time step is #pop at last time step minus (#pop at current time step minus #new births)
      D$deaths.total[2:nrow(D)] <- D$mod.pop[1:(nrow(D)-1)] - (D$mod.pop[2:nrow(D)] - D$births.total[2:nrow(D)])
      D$infectious.laststep[2:nrow(D)] <- (D$inf.sub[1:(nrow(D)-1)] + D$inf.adult[1:(nrow(D)-1)]) / (D$pos.sub[1:(nrow(D)-1)] + D$pos.adult[1:(nrow(D)-1)])
      
      D <- D %>%
        mutate(birth.rate = births.total / mod.pop,
               death.rate = deaths.total / mod.pop)
      
    } else {
      D <- D %>% 
        group_by(time) %>% 
        mutate(
          mod.pop    = sum(D_all),
          sub.pop    = sum(Dsub_all),
          adult.pop  = sum(D_adult_all + J_all),
          prev.sub   = sum(Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2) / sub.pop,
          prev.adult = sum(D_adult_h + D_adult_v + J_h + J_v) / adult.pop,
          seroprev.sub   = sum(Dsub_h1 + Dsub_h2 + Dsub_v1 + Dsub_v2 + Dsub_r1 + Dsub_r2 + Dsub_m) / sub.pop,
          seroprev.adult = sum(D_adult_h + D_adult_v + D_adult_r + J_h + J_v + J_r) / adult.pop
        )
    }
    
    # Return the model output
    return(D)
    
  })
}


# AUXILIARY FUNCTIONS ----

# Function to calculate horizontal force of infection using RK4 integration
GetFoI <- function(non_inf_size, inf_size, beta, dens_cont = FALSE) {
  
  # non_inf_size: no. susceptible individuals
  # inf_size: no. infected individuals
  # beta: transmission rate
  # dens_cont: toggle for density-dependent contact rate
  
  ## Calculate population density specific mean contact rates
  if (dens_cont) {
    # Calculate population density (using MOSA's 3 ha field site as the area)
    Nd <- (non_inf_size + inf_size) / 3
    # Calculate density dependent mean no. contacts using fitted sigmoid function (https://nsojournals.onlinelibrary.wiley.com/doi/abs/10.1111/oik.03623)
    n_cont <- 2.13 / ( 1 + exp( -0.10 * (Nd - 50.60) ) )
    scaled_n_cont <- n_cont / 2.13
  } else {
    scaled_n_cont <- 1
  }
  
  ## Calculate RK4 slopes
  
  # Time step t: Calculate rate of new infections
  k_1 <- (beta * scaled_n_cont) * non_inf_size * inf_size  # slope = initial value of derivative 
  
  # Time step t+1/2: first order estimate of slope at midpoint of t & t+1
  k_2 <- (beta * scaled_n_cont) * (non_inf_size - k_1 / 2) * (inf_size + k_1 / 2)
  k_2 <- max(k_2, 0)
  
  # Time step t+1/2: second order estimate of slope at midpoint of t & t+1
  k_3 <- (beta * scaled_n_cont) * (non_inf_size - k_2 / 2) * (inf_size + k_2 / 2)
  k_3 <- max(k_3, 0)
  
  # Time step t+1: estimate using slope of k3
  k_4 <- (beta * scaled_n_cont) * (non_inf_size - k_3) * (inf_size + k_3)
  k_4 <- max(k_4, 0)
  
  # Estimate no. horizontal transmission events btwn time t & t+1 as weighted average of slopes
  E_total <- (1/6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
  
  ## Derive & return the force of infection
  h_inf <- E_total / non_inf_size  # per susceptible individual, hence /S
  
  return(h_inf)
}


# Function to count how many times each element of values occurs in vector
dens <- function(vector, values) {
  
  # vector: A vector of data values
  # values: A vector of unique values (categories) to check frequencies for
  
  output <- rep(0, times = length(values))
  for (i in 1:length(values)) {
    output[i] <- length(which(vector == values[i])) 
  }
  
  # Returns numeric vector of counts, where each entry corresponds to the frequency 
  # of the matching element in values
  return(output)
}
