library(tbeptools)
library(tidyverse)
library(here)
library(sf)

# NAD83(HARN) / Florida West (ftUS)
# this is the projection for the seagrass segment layer from the district
prj <- 2882

# patchy, continuous, and floating algae
flcat <-c('9113', '9116', '9121')

# sg management areas
mngs <- sgmanagement %>% 
  st_transform(crs = prj) %>% 
  mutate(
    mngacre = as.numeric(units::set_units(st_area(.), 'acres'))
  )

# all zipped files on amazon s3
# downloaded from here https://data-swfwmd.opendata.arcgis.com/
sgyrs <- c('88', '90', '92', '94', '96', '99', '01', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22') 
sgyrs <- ifelse(as.numeric(sgyrs) >= 88, paste0('19', sgyrs), paste0('20', sgyrs))

# fim stations, add sgyear as interval within seagrass coverage years, right closed
# assumes coverage years are taken one year prior
fimsta <- fimstations %>% 
  mutate(
    year = str_sub(Reference, 4, 7), 
    sgyear = sgyrs[1 + findInterval(year, as.numeric(sgyrs) - 1, left.open = T)]
  ) %>% 
  select(-bay_segment) %>% 
  st_transform(crs = prj)

# unique seagrass years for fim data, for iterating
fimsgyrs <- fimsta$sgyear %>% 
  unique() 

out <- NULL
for(i in seq_along(fimsgyrs)){

  cat(i, 'of', length(fimsgyrs), '\n')
  
  fimsgyr <- fimsgyrs[i]
  
  ## import seagrass file
  
  fl <- gsub('^19|^20', '', fimsgyr) %>% 
    paste0('https://swfwmd-seagrass.s3.amazonaws.com/sg', ., '.zip')
  
  # download from s3, unzip
  tmpdir <- here('data/tmp')
  tmpzip <- here('data/tmp/tmp.zip')
  dir.create(tmpdir)
  download.file(fl, destfile = tmpzip)
  unzip(tmpzip, exdir = tmpdir)
  
  # import sg shapefile
  toimp <- list.files(tmpdir, pattern = '\\.shp$', full.names = T)
  dat_raw <- st_read(toimp, quiet = T)
  
  # delete files
  unlink(tmpdir, recursive = T)
  
  if(any(c('FLUCCS_CODE', 'FLUCCS_COD') %in% names(dat_raw)))
    names(dat_raw) <- gsub('^FLUCCS\\_COD$|^FLUCCS\\_CODE$', 'FLUCCSCODE', names(dat_raw))
  
  # filter by fluccs, intersect with sg management areas, get acreages
  # 9113 is patchy, 9116 is continuous
  dat_crp <- dat_raw %>%
    st_transform(crs = prj) %>%
    dplyr::select(FLUCCSCODE) %>% 
    filter(FLUCCSCODE %in% flcat) %>%
    st_intersection(mngs) %>% 
    summarise(
      geometry = st_union(geometry), 
      .by = c('FLUCCSCODE', 'areas', 'mngacre')
    ) %>% 
    mutate(
      acres = as.numeric(units::set_units(st_area(.), 'acres'))
    )
  
  # get fim data by sgyr
  fimtmp <- fimsta %>% 
    filter(sgyear == fimsgyr)
  
  # intersect fim stations for sgyr by dat_crp
  fimint <- fimtmp %>% 
    st_intersects(., dat_crp) %>% 
    .[1:length(.)]
  fluccsyr <- lapply(fimint, function(x) ifelse(length(x) > 0, x, NA)) %>% 
    unlist() %>% 
    dat_crp[., c('FLUCCSCODE', 'acres')] %>% 
    st_set_geometry(NULL)
  
  # add fluccs intersect to fim stations, get sg management area for each station
  fimout <- fimtmp %>%
    bind_cols(fluccsyr) %>% 
    st_intersection(mngs) %>% 
    mutate(FLUCCSCODE = as.integer(FLUCCSCODE))

  row.names(fimout) <- 1:nrow(fimout)
  
  out <- bind_rows(out, fimout)
  
}

# save as RData object
fimsgdat <- out
save(fimsgdat, file = here('data/fimsgdat.RData'))

# save as csv
fimsgdat <- fimsgdat %>% 
  st_transform(crs = 4326) %>% 
  mutate(
    lon = st_coordinates(.)[, 1], 
    lat = st_coordinates(.)[, 2]
  ) %>% 
  st_set_geometry(NULL)
write.csv(fimsgdat, file = here('data/fimsgdat.csv'), row.names = F)
