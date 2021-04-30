Theta_Logistic_Model = function(N0, r, K, sigma, theta, years, random_seed){
  
  set.seed(random_seed)
  
  N_record = N = N0
  
  for(year in seq(years)){
    epsilon = rnorm(1, mean = 0, sd = sigma)
    N = round(exp(log(N) + r * (1 - (N / K) ^ theta) + epsilon))
    N_record = c(N_record, N)
    
    if(N == 0) break
  }
  
  return(tibble(year = 0:year, N = N_record))
}

Exp_Model = function(N0, r, sigma, year0, years, random_seed){
  
  set.seed(random_seed)
  
  N_record = N = N0
  
  for(year in seq(years)){
    epsilon = rnorm(1, mean = 0, sd = sigma)
    N = round(exp(log(N) + r + epsilon))
    N_record = c(N_record, N)
    
    if(N == 0) break
  }
  
  return(tibble(year = year0 + 0:years, N = N_record))
}

# Process PVA results for quasiextinction statistics
quasiextinction_analysis = function(simulation_results, parameters){
  
  invisible(list2env(parameters, envir = environment()))
  
  if(model == 'Theta_Logistic'){
    result = 
      simulation_results %>%
      filter(N <= quasiextinction_threshold) %>%
      group_by(r, K, theta, sigma, random_seed) %>%
      slice_min(year) %>%
      ungroup %>%
      group_by(r, K, theta, sigma) %>%
      mutate(cumulative_prob = rank(year) / replicates * 100) %>%
      ungroup
  }
  
  if(model == 'Exp_Model'){
    result = 
      simulation_results %>%
      filter(N <= quasiextinction_threshold) %>%
      group_by(r, sigma, random_seed) %>%
      slice_min(year) %>%
      ungroup %>%
      group_by(r, sigma) %>%
      mutate(cumulative_prob = rank(year) / replicates * 100) %>%
      ungroup
  }
  
  
  return(result)
}

# Simulate PVA based on fitted model parameters
simulate_PVA = function(parameters){
  
  invisible(list2env(parameters, envir = environment()))
  
  if(model == 'Theta_Logistic'){
    result =
      expand_grid(
        N0 = N0, 
        r = r, 
        K = K,
        theta = theta, 
        sigma = sigma, 
        years = simulation_years, 
        random_seed = 1:replicates
      ) %>%
      mutate(model = pmap(., Theta_Logistic_Model)) %>%
      unnest(cols = c(model))
  }
  
  if(model == 'Exp_Model'){
    result =
      expand_grid(
        N0 = N0, 
        r = r, 
        sigma = sigma, 
        year0 = year0,
        years = simulation_years, 
        random_seed = 1:replicates
      ) %>%
      mutate(model = pmap(., Exp_Model)) %>%
      unnest(cols = c(model))
  }
  
  
  return(result)
} 

# Plot quasiextinction analysis
plot_extinction_risk = function(input_dtf, parameters){
  
  invisible(list2env(parameters, envir = environment()))
  
  if(!is.null(focal_parameter)){
    plot =
      input_dtf %>%
      rename(focal_parm = focal_parameter) %>%
      mutate(focal_parm = as.factor(round(focal_parm, 2))) %>%
      ggplot(aes(year, cumulative_prob, color = focal_parm, group = focal_parm)) +
      geom_line() +
      labs(color = focal_parameter)  
  }
  
  if(is.null(focal_parameter)){
    plot =
      input_dtf %>%
      ggplot(aes(year, cumulative_prob)) +
      geom_line()
  }
  
  plot = 
    plot +
    ylab('cumulative probability of quasiextinction [%]') +
    theme(aspect.ratio = 1) +
    ggtitle(paste('Quasiextinction threshold: N =', quasiextinction_threshold))
  
  return(plot)
  
}

# Wrapper
PVA_routine = function(parms){
  x = simulate_PVA(parms)
  y = quasiextinction_analysis(x, parms)
  z = plot_extinction_risk(y, parms)
  
  return(z)
}


