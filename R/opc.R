#' opcr
#'
#' An Optical Plankton Counter (OPC) is a tried and true piece of oceanographic equipment
#' designed to count and size small particles in the water. The raw data are stored in a
#' binary file format (typically with a .D00 file extension) designed by the manufacturer,
#' Focal Technologies. This package provides utilities for reading in D00 files, converting
#' them to a nice R format, and deriving and plotting relevant metrics.
#'
#' @docType package
#' @author Hansen Johnson (\email{hansen.johnson@@dal.ca})
#' @name opcr
NULL

# data --------------------------------------------------------------------

#' Processed OPC downcast
#'
#' A dataset containing a processed, trimmed OPC downcast (cast id `2019b_27`) that
#' was collected in the southern Gulf of St Lawrence in August, 2019.
#'
#' @format A nested tibble with 148 rows and 8 variables:
#' \describe{
#'   \item{scan}{the scan number, which iteratively increases with each data record}
#'   \item{timer}{the timer, which recorded seconds since the unit was powered on}
#'   \item{atten}{light attenuation}
#'   \item{depth}{instrument depth in meters}
#'   \item{flag}{quality control flag, with zero meaning good. See `opc_flag()` for the other definitions}
#'   \item{time}{the datetime since deployment, in UTC}
#'   \item{volume_filtered}{the volume of water that has passed through the OPC since the previous record, in cubic meters}
#'   \item{esd}{a nested list of the particle sizes, in Equivalent Spherical Diameter (ESD; mm) detected during this data record}
#' }
#' @source hansen.johnson@dal.ca
"opc"

# process -----------------------------------------------------------------

#' Read OPC raw data file
#'
#' OPC data are stored in a unique binary file format with *.D00 file
#' extension. This function was adapted from Mark Baumgartner's IDL routine
#' of the same name to read data from a *.D00 file and extract the id flags
#' and data fields.
#'
#' @param ifile path to .D00 file
#' @param nheader number of bits in header (default = 27)
#' @param tformat strptime format used to extract timestamp from header
#' @param tz timezone
#'
#' @return list containing `start_time` of the data file and vectors of the `id` flags and
#' corresponding `data` values
#'
#' @export
#'
#' @examples
read_focal_opc = function(ifile, nheader = 27, tformat = 'WHOI %m/%d/%y %H:%M:%S \n\n\r', tz = 'UTC'){

  # determine file size
  fs = file.info(ifile)$size

  # number of bytes in body
  nbody = fs-nheader

  # number of two-byte words
  n = nbody/2

  # open file connection
  f = file(ifile,"rb")

  # read in header
  header = readBin(f, 'raw', size=1, n=nheader, signed = F)

  # read in data
  words = readBin(f, 'int', size=1, n=nbody, signed = F)

  # close file
  close(f)

  # reform data into two byte words
  dim(words) = c(2,n)
  words = t(words)

  # bitwise operations to select id and data bits from each word
  id = bitwAnd(words[,1], 240) / 16
  data = as.integer(bitwAnd(words[,1],15)) * 256 + as.integer(words[,2])

  # identify bad records
  ii = 1
  flag = rep(1,n)
  while(ii < n){
    if(id[ii] == 11 | id[ii] == 14){
      skip = switch(as.character(id[ii]),'11'= 4, '14' = 13)
      jj = (ii + skip-1)
      if(jj > (n-1)){jj = n-1}
      flag[ii:jj] = 0
      ii = jj
    }
    ii = ii+1
  }

  # select only good records
  good = which(flag == 1)
  if(length(good)>0){
    id = id[good]
    data = data[good]
  }

  # convert header to timestamp
  txt = readBin(header, what = 'character')
  start_time = as.POSIXct(txt, format=tformat, tz=tz)

  # return data
  return(
    list(
      'start_time' = start_time,
      'id' = id,
      'data' = data
    )
  )
}

#' Convert digital size
#'
#' Apply the calibration function from the Focal OPC manual (Appendix A)
#' to convert from raw 1-4096 digital size bins recorded by the OPC to
#' Equivalent Spherical Diameter (ESD) in mm. This is called within
#' `convert_single_opc()`. It was modified from Mark Baumgartner's IDL
#' routine of the same name.
#'
#' @param ds numeric
#'
#' @return numeric
#' @export
#'
#' @examples
convert_digital_size = function(ds){

  a0 = 10879.8
  a1 = 1.923
  s = sqrt(ds)
  term = a0 / ((exp(abs(3.65 - (ds / 1000.0)))) ^ a1)
  esd = (2088.76 + term + 85.85 * s) * (1 - exp(-0.0465 * s + 0.00008629 * ds))
  esd = esd / 1000.0

  return(esd)
}

#' Convert OPC data
#'
#' Read raw OPC data and convert to nice tabular format. This was modified (heavily) from
#' from Mark Baumgartner's IDL routine of the same name.
#'
#' @param ifile path to .D00 file
#'
#' @return tibble
#' @export
#'
#' @examples
convert_single_opc = function(ifile){

  # read in D00 file and extract values
  d = read_focal_opc(ifile)
  start_time = d$start_time
  id = d$id
  data = d$data

  # calibration coefficients
  depth_sp_grav = 1.027
  flow_m = 0.130
  flow_b = 0.0370
  depth_m = 250.0
  depth_b = -250.0
  fill = 9.9692099683868690e+36
  area = 0.02 * 0.25

  # define starting variables
  first_timer = 0
  atten = 4095
  depth = fill
  flow = fill

  # pre-allocate variables
  j = which(id == 3)
  rec = tibble(scan=seq(1,length(j),1),timer=0,atten=0,depth=0,flow=0,flag=0,start=0,count=0)
  j = which(id == 1)
  esd = rep(0,length(j))

  for(ii in seq_along(id)){

    if(id[ii] == 1){ # particle record

      if(first_timer == 1){
        esd[k] = convert_digital_size(data[ii])
        if(count == 0){
          rec$start[j] = k
        }
        k = k+1
        count = count+1
      }

    } else if(id[ii] == 2){ # light attenuance record

      atten = data[ii]

    } else if(id[ii] == 3){ # time record

      if(first_timer == 0){
        j = 1
        k = 1
        first_timer = 1
      } else {
        if(count == 0){
          rec$start[j] = -1
        }
        rec$count[j] = count
        j = j+1
      }

      count = 0
      rec$timer[j] = data[ii]
      rec$atten[j] = atten
      rec$depth[j] = depth
      rec$flow[j] = flow
      atten = 4095
      depth = fill

    } else if(id[ii] == 5){ # depth record

      x = 5.0 * data[ii] / 4095.0
      p = depth_m * x + depth_b
      depth = p * 0.70307 / depth_sp_grav

    } else if(id[ii] == 8){ # flow record

      flow = flow_m * data[ii] + flow_b

    }
  }

  # add last count
  if(count == 0){
    rec$start[j] = -1
  }
  rec$count[j] = count

  # define total number of records
  nrec = j+1
  n = k

  # fix timer (resets at 4095) and convert to real time
  t = c(1,diff(rec$timer))
  t[t==-4095] = 1
  rec$time = start_time+cumsum(t)

  # calculate volume filtered
  rec$volume_filtered = abs(c(diff(rec$depth),0))*area

  # add esd to rec as nested list
  rec$esd = NA
  for(ii in 1:nrow(rec)){
    if(rec$start[ii]!=-1){
      s0 = rec$start[ii]
      s1 = s0+rec$count[ii]-1
      rec$esd[ii] = list(esd[s0:s1])
    }
  }

  # remove unused columns
  rec$flow = NULL
  rec$start = NULL
  rec$count = NULL

  return(rec)
}

#' Flag bad data
#'
#' Flag bad OPC data based on several criteria. When a record is flagged
#' the `flag` value is updated to reflect the reason. The criteria and
#' associated flag values are as follows:
#'
#' 1) extreme depths - `depth`
#' 2) non-sequential timer values - `timer`
#' 3) slow descent rates (<0.3 m/s) - `slow`
#' 4) directional reversals - `reversal`
#' 5) extreme light attenuance - `attenuance`
#'
#' @param df OPC data frame
#'
#' @return OPC tibble with updated flag column
#' @export
#'
#' @examples
opc_flag = function(df){

  # calculate time, depth, and attenuation difference
  depth_diff = c(NA,diff(df$depth))
  time_diff = c(NA,diff(df$timer))
  atten_diff = c(NA,diff(df$atten))

  # calculate instantaneous speed
  speed = depth_diff / time_diff

  # flag non-sequential timer values
  df$flag[time_diff != 1] = 'timer'

  # flag slow descent speeds
  df$flag[speed < 0.3] = 'slow'

  # flag direction reversals
  df$flag[depth_diff < 0] = 'reversal'

  # flag excessive change in attenuance
  df$flag[abs(atten_diff) > 50] = 'attenuance'

  # flag excessive attenuance
  df$flag[df$atten > 1800] = 'attenuance'

  # flag extreme depth values
  df$flag[df$depth > 1e30] = 'depth'

  # apply median filter to flag depth spikes
  mf = stats::runmed(df$depth, k = 11)
  df$flag[which(abs(mf-df$depth)>3)] = 'depth'

  return(df)
}

#' Trim OPC cast
#'
#' Shiny app to select the downcast portion of an OPC cast. It also applies
#' `opc_flag()` to flag bad values and uses `opc_plot_diagnostics()` to produce
#' a helpful diagnostic plot of the selected region.
#'
#' @param df OPC tibble
#'
#' @return OPC tibble trimmed to include selected downcast
#' @export
#'
#' @examples
opc_trim = function(df){
  # shiny app to select downcast
  ui = fluidPage(
    fluidRow(
      column(width = 12,
             helpText('Click and drag to select a region. Double click inside a selected region to zoom in, or outside to reset the plot limits.', align = "center"),
             plotOutput("full", height = 400, dblclick = "plot_dblclick",
                        brush = brushOpts(id = "plot_brush", direction = "x",resetOnNew = TRUE)
                        )
      )
    ),
    fluidRow(
      column(width = 12,
             helpText('Click `Plot` to plot OPC data in the selected region. Click `Done` to trim and save the output', align = "center")
             ),
      column(width = 3, offset = 3,
             actionButton("plot", label = 'Plot',width = '100%')
             ),
      column(width = 3,
             actionButton("done", label = 'Done',width = '100%')
             )
    ),
    fluidRow(
      column(width = 12,plotOutput("diagnostics", height = 600)
      )
    )
  )

  server = function(input, output) {

    # initiate ranges
    ranges = reactiveValues(x = NULL)

    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot_dblclick, {
      brush = input$plot_brush
      if (!is.null(brush)) {
        ranges$x = c(brush$xmin, brush$xmax)
      } else {
        ranges$x = NULL
      }
    })

    # plot full time-depth series
    output$full <- renderPlot({
      ggplot(df[df$flag!='depth',],aes(x=scan,y=depth))+
        geom_path()+
        geom_point(shape = 1)+
        scale_y_reverse()+
        coord_cartesian(xlim = ranges$x,expand = FALSE)+
        labs(x = 'Scan', y = 'Depth (m)')+
        theme_bw()
    })

    # subset
    dfs = eventReactive(input$plot,{
      brush = input$plot_brush
      if (!is.null(brush)) {
        df %>% dplyr::filter(scan >= brush$xmin & scan <= brush$xmax)
      }else{
        showNotification("Plotting all data. Select a region to trim the data!", type = 'warning')
        df
      }
    })

    # plot diagnostics
    output$diagnostics <- renderPlot({
      if(nrow(dplyr::filter(dfs(),flag==0))>10){
        opc_plot_diagnostics(df = dfs(), good_only = TRUE)
      } else {
        showNotification("Few unflagged observations detected. Plotting all data instead...", type = 'warning')
        opc_plot_diagnostics(df = dplyr::filter(dfs(),flag!='depth'), good_only = FALSE)
      }
    })

    # return trimmed output
    observeEvent(input$done, {
      if(input$plot==0){
        showNotification("Plot the data to check the trimming!", type = 'warning')
      } else {
        stopApp(dfs())
      }
    })

  }

  runApp(shinyApp(ui, server), quiet = TRUE)
}

#' Process a single OPC cast
#'
#' The processing occurs in several steps:
#' 1. read in binary OPC data in .D00 file with `read_focal_opc()`
#' 2. convert to a nice tabular format with `convert_single_opc()`
#' 3. apply various quality control flags with `opc_flag()`
#' 4. use a shiny app to interactively select the downcast with `opc_trim()`
#'
#' @param ifile path to .D00 file
#'
#' @return opc tibble
#' @export
#'
#' @examples
opc_process = function(ifile){

  # convert data file
  opc = convert_single_opc(ifile = ifile)

  # flag bad depths before downcast selection
  opc = opc_flag(opc)

  # select downcast
  opc = opc_trim(opc)

  return(opc)
}

#' Process multiple OPC casts
#'
#' Process all the OPC raw data files in `data_dir` and save processed downcasts in
#' the `output_dir`. The naming convention extracts the cast number from the D00 file
#' name, and saves each downcast according to the convention:
#'(output_dir)/(cruise)_(cast).rds
#'
#' @param cruise prefix to add to each processed data file
#' @param data_dir path to raw (D00) data files
#' @param output_dir path to processed data files
#' @param overwrite overwrite processed data? (default is FALSE)
#'
#' @return OPC tibble
#' @export
#'
#' @examples
opc_process_cruise = function(cruise,data_dir,output_dir=data_dir,overwrite=FALSE){

  # create output dir
  if(!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

  # list data files
  flist = list.files(data_dir, pattern = '^OPC(\\d{3}).D00$', full.names = TRUE)

  # process all data files
  DF = vector('list',length(flist))
  for(ii in seq_along(flist)){

    # define raw data file
    ifile = flist[ii]

    # extract sample id
    sample_id = as.numeric(substr(basename(ifile), 4, 6))

    # interim file
    tfile = paste0(output_dir, cruise, '_', sample_id, '.rds')

    message('Processing cast ', sample_id, ' from cruise ', cruise)

    if(file.exists(tfile) & !overwrite){

      # read processed data
      message('Reading processed data from: ', tfile)
      DF[[ii]] = readRDS(tfile)

    } else {

      # process data
      opc = opc_process(ifile) %>%
        mutate(
          cruise = cruise,
          sample_id = sample_id
        )

      # save
      message('Saving processed data as: ', tfile)
      saveRDS(opc,tfile)
      DF[[ii]] = opc

    }
  }

  # flatten files
  df = bind_rows(DF)

  return(df)
}

# utils -------------------------------------------------------------------

#' Get x axis limits from ggplot
#'
#' @param p ggplot object
#'
#' @return numeric, x-axis range
#' @export
#'
#' @examples
get_xlims = function(p){
  ggplot_build(p)$layout$coord$limits$x
}

#' Get y axis limits from ggplot
#'
#' @param p ggplot object
#'
#' @return numeric, y-axis range
#' @export
#'
#' @examples
get_ylims = function(p){
  ggplot_build(p)$layout$coord$limits$y
}

#' Bin numeric vector
#'
#' @param x numeric vector
#' @param d bin width
#' @param bmin minimum bin (defaults to `0`)
#' @param bmax maximum bin (defaults to `max(x,na.rm=T)`)
#' @param ... additional arguments passed to `cut()`
#'
#' @return factor, indicating bin assignment
#' @export
#'
#' @examples
bin = function(x, d, bmin = 0, bmax = max(x,na.rm = T),...){

  # define bins
  bins = seq(from = bmin, to = ceiling(bmax/d)*d, by = d)

  # apply bins
  cut(x=x, breaks = bins, include.lowest = T, right = FALSE, ...)
}

#' Relabel numeric bin
#'
#' @param x factor indicating bin assignment
#'
#' @return numeric, indicating left hand side of bin interval (`a` of `(a,b]`)
#' @export
#'
#' @examples
relabel = function(x){

  # convert to text
  tmp = as.character(x)

  # remove brackets
  tmp = gsub('\\[','',tmp)
  tmp = gsub('\\]','',tmp)
  tmp = gsub('\\(','',tmp)
  tmp = gsub('\\)','',tmp)

  # select left side
  tmp = sapply(strsplit(tmp, ','), FUN = function(x){x[1]})

  # return numeric label
  as.numeric(tmp)
}

# calculate ---------------------------------------------------------------

#' Calculate OPC abundance profile
#'
#' @param df opc tibble
#' @param dz depth bin width (meters)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return tibble
#' @export
#'
#' @examples
#'
#' data(opc)
#'
opc_abundance = function(df, dz = 2, min_size = 1, max_size = 4, good_only = TRUE){

  # reject flagged values
  if(good_only){df = dplyr::filter(df,flag==0)}

  # bin by depth
  df$zbin = bin(df$depth,d=dz)

  # volume by depth bin
  vol = df %>%
    group_by(zbin, .drop = FALSE) %>%
    dplyr::summarize(
      v = sum(volume_filtered),
      .groups = 'drop'
    )

  # bugs by depth bin
  cnt = df %>%
    unnest_longer(esd) %>%
    dplyr::filter(esd >= min_size & esd <= max_size) %>%
    group_by(zbin, .drop = FALSE) %>%
    dplyr::summarize(
      c = n(),
      .groups = 'drop'
    )

  # combine and format
  out = full_join(cnt,vol,by='zbin') %>%
    transmute(
      depth = relabel(zbin),
      count = c,
      volume = v,
      concentration = c/v
    )

  # catch infinite
  out$concentration[is.infinite(out$concentration)]=NA

  return(out)
}

#' Calculate OPC biomass
#'
#' @param df opc tibble
#' @param dz depth bin width (meters)
#' @param good_only use only good (unflagged) values
#'
#' @return tibble
#' @export
#'
#' @examples
opc_biomass = function(df,dz=2,good_only=T){

  # reject flagged values
  if(good_only){df = dplyr::filter(df,flag==0)}

  # parameters
  rho = 1

  # bin by depth
  df$zbin = bin(df$depth,d=dz)

  # volume by depth bin
  vol = df %>%
    group_by(zbin, .drop = FALSE) %>%
    dplyr::summarize(
      v = sum(volume_filtered,na.rm = T),
      .groups = 'drop'
    )

  # biomass by depth bin
  bio = df %>%
    unnest_longer(esd) %>%
    dplyr::filter(esd >= 1) %>%
    group_by(zbin, .drop = FALSE) %>%
    dplyr::summarize(
      m = sum(4/3 * pi * (esd/2)^3 * rho, na.rm = T),
      .groups = 'drop'
    )

  # combine and format
  out = full_join(bio,vol,by='zbin') %>%
    transmute(
      depth = relabel(zbin),
      mass = m,
      volume = v,
      concentration = m/v
    )

  # catch infinite
  out$concentration[is.infinite(out$concentration)]=NA

  return(out)
}

#' Calculate OPC size-frequency histogram
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return tibble
#' @export
#'
#' @examples
opc_histogram = function(df,ds=0.05,min_size=1,max_size=4,good_only=T){

  # reject flagged values
  if(good_only){df = dplyr::filter(df,flag==0)}

  # get sampled volume
  volume = sum(df$volume_filtered)

  # extract particle sizes
  df = df %>%
    unnest_longer(esd) %>%
    dplyr::filter(esd >= min_size & esd <= max_size)

  # bin by size
  df$sbin = bin(df$esd,d=ds,bmin=min_size,bmax=max_size)

  # compute count per bin
  out = df %>%
    group_by(sbin,.drop = FALSE) %>%
    dplyr::summarize(
      count = n(),
      .groups = 'drop'
    ) %>%
    transmute(
      size = relabel(sbin),
      count,
      concentration = count/volume
    )

  # catch infinite
  out$concentration[is.infinite(out$concentration)]=NA

  return(out)
}

#' Calculate OPC abundance in size and depth bin matrix
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param dz depth bin width (m)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return tibble
#' @export
#'
#' @examples
opc_image = function(df, ds=0.05, dz=2, min_size = 1, max_size = 4, good_only = T){

  # reject flagged values
  if(good_only){df = dplyr::filter(df,flag==0)}

  # assign depth bins
  df$zbin = bin(df$depth,d=dz)

  # volume by depth bin
  vol = df %>%
    group_by(zbin, .drop = FALSE) %>%
    dplyr::summarize(
      volume = sum(volume_filtered),
      .groups = 'drop'
    )

  # extract sizes
  bg = df %>%
    unnest_longer(esd) %>%
    dplyr::filter(esd >= min_size & esd <= max_size)

  # assign size bins
  bg$sbin = bin(bg$esd,d=ds,bmin = min_size,bmax = max_size)

  # count in size and depth bins
  bg = bg %>%
    group_by(zbin,sbin, .drop = FALSE) %>%
    dplyr::summarize(
      count = n(),
      .groups = 'drop'
    )

  # combine and format
  out = full_join(bg,vol,by='zbin') %>%
    transmute(
      depth = relabel(zbin),
      size = relabel(sbin),
      count,
      volume,
      concentration = count/volume
    )

  # catch infinite
  out$concentration[is.infinite(out$concentration)]=NA

  return(out)
}

# plot --------------------------------------------------------------------

#' Plot OPC scan versus depth series
#'
#' @param df opc tibble
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_depth = function(df, good_only=T){

  if(good_only){df = dplyr::filter(df,flag==0)}

  ggplot(df,aes(x=scan,y=depth))+
    geom_path()+
    geom_point(shape=1)+
    scale_y_reverse(limits = c(NA,0))+
    labs(x='Scan',y='Depth (m)')+
    theme_bw()
}

#' Plot OPC flags
#'
#' @param df opc tibble
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_flags = function(df){
  good = dplyr::filter(df,flag==0)
  bad = dplyr::filter(df,flag!=0)
  txt = paste0(round(nrow(bad)/nrow(df)*100), '% points flagged (',nrow(bad),'/',nrow(df),')')
  ggplot()+
    geom_segment(data=bad,aes(x=scan+30,xend=scan+10,y=depth,yend=depth,color=flag),
                 alpha=0.7, arrow = arrow(length = unit(4,'points')))+
    geom_path(data=good,aes(x=scan,y=depth),alpha=0.7)+
    geom_point(data=good,aes(x=scan,y=depth),shape=1,alpha=0.7)+
    scale_color_manual(
      values=c('slow'='red','reversal'='blue','depth'='black','attenuance'='purple','timer'='orange'))+
    scale_y_reverse(limits = c(NA,0))+
    labs(x='Scan',y='Depth (m)', color = 'Flag',caption = txt)+
    theme_bw()+
    theme(legend.position = c(1,1),
          legend.justification = c(1,1),
          legend.background = element_rect(color='black'))
}

#' Plot OPC attenuance versus depth
#'
#' @param df opc tibble
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_attenuance = function(df, good_only=T){

  if(good_only){df = dplyr::filter(df,flag==0)}

  ggplot(df,aes(x=atten,y=depth))+
    geom_path(alpha=0.7)+
    geom_point(shape=1,alpha=0.7)+
    scale_y_reverse(limits = c(NA,0))+
    labs(x='Attenuance',y='Depth (m)')+
    theme_bw()
}

#' Compute and plot OPC abundance vs depth
#'
#' @param df opc tibble
#' @param dz depth bin width (m)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_abundance = function(df,dz=4,min_size=1,max_size=4,good_only=T){

  # count by depth bin
  d = opc_abundance(df = df, dz = dz, min_size = min_size,max_size = max_size, good_only = good_only) %>%
    dplyr::filter(!is.na(concentration))

  # construct x label
  xlab = paste0('Abundance (ind/m3)\n(',min_size,'-',max_size,' mm ESD)')

  # plot
  ggplot()+
    geom_rect(data=d,aes(xmin=0,xmax=concentration,ymin=depth,ymax=depth+dz),
              fill = 'grey', color = 'black', size = 0.2)+
    scale_y_reverse()+
    labs(y = 'Depth (m)', x = xlab)+
    theme_bw()
}

#' Compute and plot OPC biomass vs depth
#'
#' @param df opc tibble
#' @param dz depth bin width (m)
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_biomass = function(df,dz=4,good_only=T){

  # count by depth bin
  d = opc_biomass(df = df, dz = dz, good_only = good_only) %>%
    dplyr::filter(!is.na(concentration))

  # plot
  ggplot()+
    geom_rect(data=d,aes(xmin=0,xmax=concentration,ymin=depth,ymax=depth+dz),
              fill = 'grey', color = 'black', size = 0.2)+
    scale_y_reverse()+
    labs(y = 'Depth (m)', x = 'Biomass (ug/m3)')+
    theme_bw()
}

#' Compute and plot OPC size vs abundance
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_histogram = function(df,ds=0.05,min_size=1,max_size=4,good_only=T){

  # compute histogram
  d = opc_histogram(df,ds=ds,min_size=min_size,max_size=max_size,good_only=good_only) %>%
    dplyr::filter(!is.na(concentration))

  # plot
  ggplot()+
    geom_rect(data=d,aes(xmin=size,xmax=size+ds,ymin=0,ymax=count),
              fill = 'grey', color = 'black', size = 0.2)+
    coord_cartesian(xlim = c(min_size,max_size+ds))+
    labs(x = 'Equivalent Spherical Diameter (mm)', y = 'Abundance (ind/m2)')+
    theme_bw()
}

#' Compute and plot image-style plot of OPC abundance in size and depth bins
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param dz depth bin width (m)
#' @param min_size minimum particle ESD (mm)
#' @param max_size maximum particle ESD (mm)
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_image = function(df, ds=0.05, dz=2, min_size = 0, max_size = 4, good_only = T){

  # compute image
  d = opc_image(df=df, ds=ds, dz=dz, min_size = min_size, max_size = max_size, good_only = good_only)

  # plot image
  ggplot()+
    geom_rect(data=dplyr::filter(d,concentration>0),aes(xmin=size,xmax=size+ds,ymin=depth,ymax=depth+dz,fill=concentration))+
    scale_fill_viridis_c(guide = guide_colourbar(barwidth = 15,title.position = 'top',title.hjust = 0.5))+
    scale_y_reverse()+
    labs(x = 'Equivalent Spherical Diameter (mm)', y = 'Depth (m)', fill = 'Abundance (ind/m3)')+
    coord_cartesian(xlim=c(min_size,max_size),ylim=c(max(d$depth,na.rm = T)+dz,0))+
    theme_bw()+
    theme(legend.position = 'bottom')
}

#' Plot OPC histogram, image, and profile with shared axes
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param dz depth bin width (m)
#' @param min_size minimum particle ESD (mm) for full plot
#' @param max_size maximum particle ESD (mm) for full plot
#' @param amin_size minimum particle ESD (mm) for abundance profile
#' @param amax_size maximum particle ESD (mm) for abundance profile
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_multipanel = function(df, ds = 0.05,dz = 2, min_size = 0, max_size = 4,
                               amin_size = 1, amax_size = 2, good_only = T){

  # make plots
  p_hist = opc_plot_histogram(df = df, ds = ds, min_size = min_size, max_size = max_size, good_only = good_only)
  p_abund = opc_plot_abundance(df = df, dz = dz, min_size = amin_size, max_size = amax_size, good_only = good_only)
  p_img = opc_plot_image(df = df, ds = ds, dz = dz, min_size = min_size, max_size = max_size, good_only = good_only)

  # extract plot limits
  ylims = abs(get_ylims(p_img))
  xlims = get_xlims(p_img)

  # align and format
  suppressMessages({
    p_hist_f = p_hist + coord_cartesian(xlim = xlims)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
    p_abund_f = p_abund + coord_cartesian(ylim = ylims)+
      theme(panel.border = element_rect(fill=NA,color='blue'),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())
    p_img_f = p_img +
      geom_rect(aes(xmin=amin_size,xmax=amax_size,ymin=ylims[1],ymax=ylims[2]),fill=NA,color='blue')
  })

  # combine
  wrap_plots(p_hist_f, plot_spacer(), p_img_f, p_abund_f, widths = c(3, 1), heights = c(1, 2))
}

#' Plot a suite of OPC diagnostics
#'
#' @param df opc tibble
#' @param ds particle size bin width (mm)
#' @param dz depth bin width (m)
#' @param min_size minimum particle ESD (mm) for full plot
#' @param max_size maximum particle ESD (mm) for full plot
#' @param amin_size minimum particle ESD (mm) for abundance profile
#' @param amax_size maximum particle ESD (mm) for abundance profile
#' @param good_only use only good (unflagged) values
#'
#' @return ggplot
#' @export
#'
#' @examples
opc_plot_diagnostics = function(df, dz = 2, ds = 0.05, min_size = 0, max_size = 4,
                                amin_size = 1, amax_size = 2, good_only = T){

  # make plots
  p_flag = opc_plot_flags(df=df)
  p_atten = opc_plot_attenuance(df=df,good_only = good_only)
  p_hist = opc_plot_histogram(df=df, ds = ds, min_size = min_size, max_size = max_size, good_only = good_only)
  p_abund = opc_plot_abundance(df=df, dz = dz, min_size = amin_size, max_size = amax_size, good_only = good_only)
  p_img = opc_plot_image(df=df, dz = dz, min_size = min_size, max_size = max_size, good_only = good_only)

  # extract plot limits
  ylims = abs(get_ylims(p_img))
  xlims = get_xlims(p_img)

  # align and format
  suppressMessages({
    p_flag_f = p_flag + coord_cartesian(ylim = ylims)
    p_atten_f = p_atten + coord_cartesian(ylim = ylims)+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())
    p_img_f = p_img +
      geom_rect(aes(xmin=amin_size,xmax=amax_size,ymin=ylims[1],ymax=ylims[2]),fill=NA,color='blue')+
      theme(legend.position = 'none',
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank())
    p_hist_f = p_hist + coord_cartesian(xlim = xlims)+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
    p_abund_f = p_abund +
      scale_y_reverse(position='right')+
      coord_cartesian(ylim = ylims)+
      theme(panel.border = element_rect(fill=NA,color='blue'))
  })

  # combine
  wrap_plots(
    plot_spacer(),plot_spacer(),p_hist_f,plot_spacer(),
    p_flag_f, p_atten_f, p_img_f, p_abund_f,
    nrow=2,ncol=4,
    widths = c(2,1,2,1), heights = c(1, 2))
}
