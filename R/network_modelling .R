# library(tidyverse)
# library(statnet)
# library(tergm)
# library(ergm.multi)
# library(network)
# library(patchwork)
# library(modelsummary)
# library(geosphere)
# library(sna)

###############################################################

load_TIST_data <- function(file_path){
  if(!file.exists(file_path)){
    stop("File not found: ", file_path)
  }
  
  TistDat <- read.csv(file_path, stringsAsFactors = FALSE) |>
    tidyr::drop_na() |>
    dplyr::distinct()
  
  return(TistDat)
}

##################################################################################################
################################################################################################

########################################################################################################
build_site_networks <- function(tist_data,
                                site,
                                distance_threshold = 2000) {
  ############################################################
  # Filter site
  ############################################################
  
  site_data <- tist_data %>%
    filter(Proj_Area == site)
  
  
  if(nrow(site_data) == 0){
    stop("No records found for site: ", site)
  }
  
  
  ############################################################
  # Sets
  ############################################################
  
  all_farmers   <- sort(unique(site_data$Farmer_ID))
  all_groups    <- sort(unique(site_data$Group_ID))
  all_locations <- sort(unique(site_data$Location_ID))
  
  
  ############################################################
  # Attributes
  ############################################################
  
  group_attr <- site_data %>%
    select(Group_ID, Years_since_reg, Cluster_ID) %>%
    distinct() %>%
    mutate(
      Years_since_reg = as.numeric(Years_since_reg)
    )
  
  
  location_attr <- site_data %>%
    group_by(Location_ID) %>%
    summarise(
      Admin_Districts = first(Admin_Districts),
      Subcounty = first(Subcounty),
      Exposure = mean(Exposure, na.rm = TRUE),
      Dist_To_Forest_m = mean(Dist_To_Forest_m, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ############################################################
  # 1. FARMER-GROUP NETWORK
  ############################################################
  
  FG_df <- site_data %>%
    mutate(Count = 1) %>%
    pivot_wider(
      id_cols = Farmer_ID,
      names_from = Group_ID,
      values_from = Count,
      values_fn = sum,
      values_fill = 0
    ) %>%
    column_to_rownames("Farmer_ID") %>%
    as.matrix()
  
  
  # Add missing farmers
  
  missing_farmers <- setdiff(
    all_farmers,
    rownames(FG_df)
  )
  
  if(length(missing_farmers) > 0){
    
    FG_df <- rbind(
      FG_df,
      matrix(
        0,
        nrow = length(missing_farmers),
        ncol = ncol(FG_df),
        dimnames = list(
          missing_farmers,
          colnames(FG_df)
        )
      )
    )
  }
  
  
  # Add missing groups
  
  missing_groups <- setdiff(
    all_groups,
    colnames(FG_df)
  )
  
  
  if(length(missing_groups) > 0){
    
    FG_df <- cbind(
      FG_df,
      matrix(
        0,
        nrow = nrow(FG_df),
        ncol = length(missing_groups),
        dimnames = list(
          rownames(FG_df),
          missing_groups
        )
      )
    )
  }
  
  
  FG_df <- FG_df[
    all_farmers,
    all_groups
  ]
  
  
  FG <- network(
    FG_df,
    matrix.type = "bipartite",
    bipartite = nrow(FG_df),
    directed = FALSE
  )
  
  
  vertex_names_FG <- network.vertex.names(FG)
  
  
  for(i in seq_along(vertex_names_FG)){
    
    vname <- vertex_names_FG[i]
    
    
    if(vname %in% group_attr$Group_ID){
      
      row <- group_attr %>%
        filter(Group_ID == vname)
      
      set.vertex.attribute(
        FG,
        "type",
        "Group",
        v=i
      )
      
      set.vertex.attribute(
        FG,
        "Years_since_reg",
        row$Years_since_reg[1],
        v=i
      )
      
      set.vertex.attribute(
        FG,
        "Cluster",
        row$Cluster_ID[1],
        v=i
      )
      
    } else {
      
      set.vertex.attribute(
        FG,
        "type",
        "Farmer",
        v=i
      )
    }
  }
  
  
  
  ############################################################
  # 2. FARMER-LOCATION NETWORK
  ############################################################
  
  FL_df <- site_data %>%
    group_by(Farmer_ID, Location_ID) %>%
    filter(Trees >= 1) %>%
    summarise(
      Total_Trees = sum(Trees),
      .groups="drop"
    ) %>%
    pivot_wider(
      id_cols = Farmer_ID,
      names_from = Location_ID,
      values_from = Total_Trees,
      values_fill = 0
    ) %>%
    column_to_rownames("Farmer_ID") %>%
    as.matrix()
  
  
  # Add missing farmers
  
  missing_farmers <- setdiff(
    all_farmers,
    rownames(FL_df)
  )
  
  if(length(missing_farmers)>0){
    
    FL_df <- rbind(
      FL_df,
      matrix(
        0,
        nrow=length(missing_farmers),
        ncol=ncol(FL_df),
        dimnames=list(
          missing_farmers,
          colnames(FL_df)
        )
      )
    )
  }
  
  
  # Add missing locations
  
  missing_locations <- setdiff(
    all_locations,
    colnames(FL_df)
  )
  
  
  if(length(missing_locations)>0){
    
    FL_df <- cbind(
      FL_df,
      matrix(
        0,
        nrow=nrow(FL_df),
        ncol=length(missing_locations),
        dimnames=list(
          rownames(FL_df),
          missing_locations
        )
      )
    )
  }
  
  
  FL_df <- FL_df[
    all_farmers,
    all_locations
  ]
  
  
  FL <- network(
    FL_df,
    matrix.type="bipartite",
    bipartite=nrow(FL_df),
    directed=FALSE
  )
  
  
  vertex_names_FL <- network.vertex.names(FL)
  
  
  for(i in seq_along(vertex_names_FL)){
    
    vname <- vertex_names_FL[i]
    
    
    if(vname %in% location_attr$Location_ID){
      
      row <- location_attr %>%
        filter(Location_ID == vname)
      
      set.vertex.attribute(
        FL,
        "type",
        "Location",
        v=i
      )
      
      set.vertex.attribute(
        FL,
        "Admin_Districts",
        row$Admin_Districts[1],
        v=i
      )
      
      set.vertex.attribute(
        FL,
        "Subcounty",
        row$Subcounty[1],
        v=i
      )
      
      set.vertex.attribute(
        FL,
        "Exposure",
        row$Exposure[1],
        v=i
      )
      
      set.vertex.attribute(
        FL,
        "Dist_To_Forest_m",
        row$Dist_To_Forest_m[1],
        v=i
      )
      
    } else {
      
      set.vertex.attribute(
        FL,
        "type",
        "Farmer",
        v=i
      )
    }
  }
  
  
  
  ############################################################
  # 3. LOCATION-LOCATION NETWORK
  ############################################################
  
  location_coords <- site_data %>%
    select(
      Location_ID,
      longitude,
      latitude
    ) %>%
    distinct()
  
  
  coords_matrix <- as.matrix(
    location_coords[,c("longitude","latitude")]
  )
  
  
  dist_matrix <- distm(
    coords_matrix,
    fun = distHaversine
  )
  
  
  LL_df <- ifelse(
    dist_matrix <= distance_threshold,
    1,
    0
  )
  
  
  diag(LL_df) <- 0
  
  
  rownames(LL_df) <- location_coords$Location_ID
  colnames(LL_df) <- location_coords$Location_ID
  
  
  LL <- network(
    LL_df,
    directed=FALSE,
    matrix.type="adjacency"
  )
  
  
  
  ############################################################
  # 4. MULTILAYER NETWORK
  ############################################################
  
  FGL <- Layer(
    fg = FG,
    fl = FL
  )
  
  
  
  ############################################################
  # Replace missing attributes
  ############################################################
  
  for(attr in c(
    "Subcounty",
    "Years_since_reg",
    "Admin_Districts",
    "Dist_To_Forest_m",
    "Exposure",
    "Cluster",
    "Group_ID"
  )){
    
    vals <- FGL %v% attr
    
    
    if(is.numeric(vals)){
      
      vals[is.na(vals)] <- 0
      
    } else {
      
      vals[is.na(vals)] <- "Unknown"
      
    }
    
    
    set.vertex.attribute(
      FGL,
      attr,
      vals
    )
  }
  
  
  
  ############################################################
  # Return
  ############################################################
  
  return(
    list(
      
      site = site,
      distance_threshold = distance_threshold,
      
      FG = FG,
      FL = FL,
      LL = LL,
      FGL = FGL,
      
      FG_df = FG_df,
      FL_df = FL_df,
      LL_df = LL_df,
      
      group_attr = group_attr,
      location_attr = location_attr
      
    )
  )
}

##################################################################################################################

fit_site_model <- function(site_networks,
                           model_formula = NULL,
                           control = control.ergm()) {
  
  
  ############################################################
  # Extract networks
  ############################################################
  
  FGL <- site_networks$FGL
  LL  <- site_networks$LL
  
  
  ############################################################
  # Default model
  ############################################################
  
  if(is.null(model_formula)){
    
    model_formula <- FGL ~
      L(~edges, ~fg) +
      L(~edges, ~fl) +
      L(~b2cov("Years_since_reg"), ~fg) +
      L(~b2cov("Exposure"), ~fl) +
      L(~b2cov("Dist_To_Forest_m"), ~fl) +
      L(~Project(~edgecov(LL),2),~fl)
    
  }
  
  
  message(
    "Fitting ERGM: ",
    site_networks$site,
    " | distance = ",
    site_networks$distance_threshold,
    " m"
  )
  
  
  model <- ergm(
    model_formula,
    control = control
  )
  
  
  return(
    list(
      model = model,
      site = site_networks$site,
      distance_threshold = site_networks$distance_threshold,
      formula = model_formula,
      networks = site_networks
    )
  )
}

################################################################################################

evaluate_overlap <- function(site_model,
                             nsim = 1000,
                             seed = 123) {
  
  ################################################################################
  # Extract objects
  ################################################################################
  
  fitted_model <- site_model$model
  
  site_networks <- site_model$networks
  
  FG_df <- site_networks$FG_df
  FL_df <- site_networks$FL_df
  LL_df <- site_networks$LL_df
  
  ################################################################################
  # Observed H6 statistics
  ################################################################################
  
  FF_group <- FG_df %*% t(FG_df)
  diag(FF_group) <- 0
  
  FF_spatial <- FL_df %*% LL_df %*% t(FL_df)
  diag(FF_spatial) <- 0
  
  ## Sanity check
  stopifnot(
    identical(
      rownames(FF_group),
      rownames(FF_spatial)
    )
  )
  
  Obs_gcor <- gcor(
    FF_group,
    FF_spatial
  )
  
  Obs_overlap <- sum(
    (FF_group != 0) &
      (FF_spatial != 0)
  )
  
  cat(
    "Observed gcor:",
    Obs_gcor,
    "\n"
  )
  
  cat(
    "Observed overlap:",
    Obs_overlap,
    "\n"
  )
  
  ################################################################################
  # Simulation
  ################################################################################
  
  set.seed(seed)
  
  sim_nets <- simulate(
    fitted_model,
    nsim = nsim,
    output = "network"
  )
  
  sim_gcor <- numeric(nsim)
  sim_overlap <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    
    sim_net <- sim_nets[[i]]
    
    A <- uncombine_network(sim_net)
    
    FL_sim <- A$fl
    FG_sim <- A$fg
    
    fl_type <- FL_sim %v% "type"
    fg_type <- FG_sim %v% "type"
    
    farmers_fl <- which(fl_type == "Farmer")
    locations  <- which(fl_type == "Location")
    
    farmers_fg <- which(fg_type == "Farmer")
    groups     <- which(fg_type == "Group")
    
    FL_msim <- as.matrix.network(FL_sim)[
      farmers_fl,
      locations,
      drop = FALSE
    ]
    
    FG_msim <- as.matrix.network(FG_sim)[
      farmers_fg,
      groups,
      drop = FALSE
    ]
    
    ####################################################################
    # Ensure farmer ordering is identical
    ####################################################################
    
    fl_farmer_names <- network.vertex.names(FL_sim)[farmers_fl]
    fg_farmer_names <- network.vertex.names(FG_sim)[farmers_fg]
    
    if (!identical(
      fl_farmer_names,
      fg_farmer_names
    )) {
      
      common <- intersect(
        fl_farmer_names,
        fg_farmer_names
      )
      
      FL_msim <- FL_msim[
        match(common, fl_farmer_names),
        ,
        drop = FALSE
      ]
      
      FG_msim <- FG_msim[
        match(common, fg_farmer_names),
        ,
        drop = FALSE
      ]
      
    }
    
    ####################################################################
    # Project farmer networks
    ####################################################################
    
    FFg_sim <- FG_msim %*% t(FG_msim)
    diag(FFg_sim) <- 0
    
    FFs_sim <- FL_msim %*% LL_df %*% t(FL_msim)
    diag(FFs_sim) <- 0
    
    ####################################################################
    # Statistics
    ####################################################################
    
    sim_gcor[i] <- gcor(
      FFg_sim,
      FFs_sim
    )
    
    sim_overlap[i] <- sum(
      (FFg_sim > 0) &
        (FFs_sim > 0)
    )
    
  }
  
  ################################################################################
  # CORRECTED p-values
  ################################################################################
  
  p_gcor <- mean(
    sim_gcor >= Obs_gcor
  )
  
  p_overlap <- mean(
    sim_overlap >= Obs_overlap
  )
  
  ################################################################################
  # Return
  ################################################################################
  
  return(
    
    list(
      
      site = site_networks$site,
      
      distance_threshold =
        site_networks$distance_threshold,
      
      observed = list(
        
        gcor = Obs_gcor,
        
        overlap = Obs_overlap
        
      ),
      
      simulated = list(
        
        gcor = sim_gcor,
        
        overlap = sim_overlap
        
      ),
      
      p.values = list(
        
        gcor = p_gcor,
        
        overlap = p_overlap
        
      ),
      
      nsim = nsim,
      
      seed = seed
      
    )
    
  )
  
}
###############################################################################################
summarise_model_results <- function(model_obj,
                                    simulation_obj){
  
  fit <- model_obj$model
  
  list(
    
    coefficients =
      summary(fit)$coefficients,
    
    model_statistics =
      tibble(
        AIC = AIC(fit),
        BIC = BIC(fit),
        LogLik = as.numeric(logLik(fit))
      ),
    
    overlap =
      tibble(
        observed_gcor =
          simulation_obj$observed$gcor,
        
        simulated_gcor =
          mean(simulation_obj$simulated$gcor),
        
        p_gcor =
          simulation_obj$p.values$gcor,
        
        observed_overlap =
          simulation_obj$observed$overlap,
        
        simulated_overlap =
          mean(simulation_obj$simulated$overlap),
        
        p_overlap =
          simulation_obj$p.values$overlap
      )
    
  )
  
}