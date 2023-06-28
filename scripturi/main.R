library(tidyverse)
library(furrr)
plan(multisession, workers=4)
theme_set(hrbrthemes::theme_ipsum_tw())


iciot <- map(1995:2018, ~read_csv(citeste(.x)) %>% rename(prod=...1))

make_A <- function(icio, selection="ROU", what="A"){
  icio <- icio %>% filter(!grepl("TAXSUB|VALU", prod)) %>%
          select(-TOTAL) %>% select(-matches("HFCE|NPISH|GGFC|GFCF|INVNT|DPABR"))
  output <- icio %>% filter(prod=="OUTPUT") %>% select(-1) %>% unlist()
  output <- diag(1/output)
  output[is.infinite(output)] <- 0
  output[is.nan(output)] <- 0
  Z <- icio %>% filter(prod!="OUTPUT") %>% column_to_rownames("prod") %>% as.matrix()
  colnames(output) <- colnames(Z)
  rownames(output) <- rownames(Z)
  if(what=='Z'){
    return(Z)
  }
  if(selection=='all'){
  A <- Z %*% output
  } else {
    A <- Z[grepl(selection, rownames(Z)), 
           grepl(selection, colnames(Z))] %*% output[grepl(selection, rownames(output)), 
                                                     grepl(selection, colnames(output))]
  }
return(A)
}

make_L <- function(A){
  L <- solve((diag(1, nrow(A)) - A), diag(nrow(A)))
  colnames(L) <- colnames(A)
  rownames(L) <- rownames(A)
  return(L)
}

make_M <- function(icio, selection="ROU"){
  fdc <- "HFCE|NPISH|GGFC|GFCF|INVNT|DPABR"
  
  icio <- icio %>% filter(!grepl("TAXSUB|VALU", prod)) %>%
    select(-TOTAL) %>% select(-matches(fdc))
  
  output <- icio %>% filter(prod=="OUTPUT") %>% select(-1) %>% unlist()
  output <- diag(1/output)
  output[is.infinite(output)] <- 0
  output[is.nan(output)] <- 0
  
  Z <- icio %>% filter(prod!="OUTPUT") %>% 
       column_to_rownames("prod") %>% as.matrix()
  colnames(output) <- colnames(Z)
  rownames(output) <- rownames(Z)
  
  A <- Z[!grepl(selection, rownames(Z)), 
         grepl(selection, colnames(Z))] %*% output[grepl(selection, rownames(output)), 
                                                     grepl(selection, colnames(output))]
return(A)
}


make_V <- function(icio, selection="ROU"){
  fdc <- "HFCE|NPISH|GGFC|GFCF|INVNT|DPABR"
  
  icio <- icio %>% filter(!grepl("TAXSUB", prod)) %>%
    select(-TOTAL) %>% select(-matches(fdc))
  
  output <- icio %>% filter(prod%in%c("VALU", "OUTPUT")) %>% 
            select(prod, matches(selection)) %>%
            column_to_rownames("prod") %>% as.matrix()
  VA <- diag(output[1,]/output[2,])
  VA[is.infinite(VA)] <- 0
  VA[is.nan(VA)] <- 0
  colnames(VA) <- colnames(output)
  rownames(VA) <- colnames(VA)
  return(VA)
}

make_f <- function(icio, selection="ROU"){
  index <- paste0("^", selection, "_", "[0-9]{2}")
  f <- icio %>% filter(!grepl("TAXSUB|VALU|OUTPUT", prod)) %>% 
       select(-TOTAL) %>% select(-matches(index)) %>%
       column_to_rownames("prod") %>% as.matrix()  
  
  f[grepl(selection, rownames(f)), ]
}

make_I <- function(icio, selection="ROU"){
  f <- icio %>% filter(!grepl("TAXSUB|VALU|OUTPUT", prod)) %>% 
    select(prod, matches("HFCE|NPISH|GGFC|GFCF|INVNT|DPABR")) %>%
    column_to_rownames("prod") %>% as.matrix()  
  f[!grepl(selection, rownames(f)), grepl(selection, colnames(f))]
}

agrega <- function(matrice, axa=1, what, selection='ROU'){
  if(axa==1){
   if(what=='sectoare'){
      x <- as_tibble(matrice, rownames='prod') %>% 
        gather(2:ncol(.), key="use", value='valoare') %>% 
        group_by(use, prod = gsub("^[A-Z0-9]{3}\\_", "", prod)) %>% 
        summarise(valoare = sum(valoare, na.rm=TRUE)) %>% ungroup() %>% 
        spread(use, valoare) %>% column_to_rownames("prod") %>% 
        as.matrix()
      colnames(x) <- gsub("^.+\\_", "INT_", colnames(x))
      
    } else if(what=='tari'){
      x <- as_tibble(matrice, rownames='prod') %>% 
        gather(2:ncol(.), key="use", value='valoare') %>% 
        group_by(use, prod = gsub("\\_[0-9T]+$", "", prod)) %>% 
        summarise(valoare = sum(valoare, na.rm=TRUE)) %>% ungroup() %>% 
        spread(use, valoare) %>% column_to_rownames("prod") %>% 
        as.matrix()
      
    } else {
      stop("Alegeţi sectoare sau ţări !")
    }
  } else if(axa==2){
    if(what=='categorii'){
      x <- as_tibble(matrice, rownames='prod') %>% 
        gather(2:ncol(.), key="use", value='valoare') %>% 
        group_by(prod, use = if_else(grepl(selection, use), 
                                     use, "Exporturi")) %>% 
        summarise(valoare = sum(valoare, na.rm=TRUE)) %>% ungroup() %>% 
        spread(use, valoare) %>% column_to_rownames("prod") %>% 
        as.matrix()}
    else if(what=='sectoare'){
      x <- as_tibble(matrice, rownames='prod') %>% 
        gather(2:ncol(.), key="use", value='valoare') %>% 
        group_by(prod, use = if_else(grepl(selection, use), use, gsub("^[A-Z0-9]{3}\\_", "", use))) %>% 
        summarise(valoare = sum(valoare, na.rm=TRUE)) %>% ungroup() %>% 
        spread(use, valoare) %>% column_to_rownames("prod") %>% 
        as.matrix()
      colnames(x) <- gsub("^.+\\_", "INT_", colnames(x))
      
      
    } else if(what=='tari'){
      x <- as_tibble(matrice, rownames='prod') %>% 
        gather(2:ncol(.), key="use", value='valoare') %>% 
        group_by(prod, use = gsub("\\_[0-9T]+$", "", use)) %>% 
        summarise(valoare = sum(valoare, na.rm=TRUE)) %>% ungroup() %>% 
        spread(use, valoare) %>% column_to_rownames("prod") %>% 
        as.matrix()
      
    } else {
      stop("Alegeţi sectoare sau ţări !")
    }
  } else { stop("Alegeţi axa 1 sau 2!")
  }
  return(x)
}

aduna <- function(I, M){
  as_tibble(I, rownames='prod') %>% 
    gather(2:ncol(.), key='sectoare', value='valoare') %>%
    bind_rows(as_tibble(M, rownames="prod") %>% 
                gather(2:ncol(.), key='categorii', value='valoare')) %>%
    group_by(prod, categorii) %>% summarise(valoare=sum(valoare, na.rm=TRUE)) %>%
    ungroup() %>% spread(categorii, valoare) %>% column_to_rownames("prod") %>%
    as.matrix()
}

#de incredere
importuri <- function(icio, selection){
  i <- make_M(icio, selection = selection)%*%make_L(make_A(icio, selection=selection))%*%make_f(icio, selection=selection)
  dir <- make_I(icio, selection=selection)
  index <- colnames(i)%in%colnames(dir)
  mat <- matrix(data=0, ncol=length(index), nrow=nrow(i))
  mat[,index] <- dir
  rez <- i + mat
  return(rez)
}

eastern <- list()
for(i in c("ROU", "CZE", "HUN", "BGR", "POL", "SVK")){
eastern[[i]] <- future_map(iciot, ~importuri(icio=.x, selection=i) %>%
                       agrega2(., by='categorii', selection=i) %>%
                       agrega(., axa=1, what='sectoare', selection=i))

eastern[[i]] <- map2_dfr(eastern[[i]],1995:2018, ~.x %>% as_tibble(rownames='sectoare') %>% 
                                         mutate(year=.y))

eastern[[i]] <- eastern[[i]] %>% mutate(Exporturi=Exporturi+INT_DPABR,
                              Gospodării = INT_HFCE+INT_NPISH,
                              Stat = INT_GGFC,
                              Investiţii = INT_GFCF) %>%
           select(sectoare, Exporturi, Gospodării, Stat, Investiţii, year) %>%
           gather(Exporturi:Investiţii, key='component', value='valoare')
}

eastern <- map2_dfr(eastern, c("ROU", "CZE", "HUN", "BGR", "POL", "SVK"), 
                   ~.x %>% mutate(country=.y))



aranjeaza <- function(icio, selection){
i <- make_M(icio, selection = selection)%*%make_L(make_A(icio, selection=selection))%*%make_f(icio, selection=selection)
i <- agrega(i, axa=2, what='sectoare', selection=selection)

rez <- aduna(I=make_I(icio, selection=selection), M=i) %>% 
  agrega(., axa=1, what='sectoare', selection=selection) %>% 
  as_tibble(rownames='prod') 
return(rez)
}


importuri_tara <- function(country){
rez <- future_map(c(1:24),~aranjeaza(icio=iciot[[.x]], selection=country) %>% 
                    mutate(year=.x+1994)) %>% bind_rows() %>%
                    group_by(year) %>% summarise(across(2:57, sum)) %>%
                    ungroup() %>% gather(2:ncol(.), key='name', value='imp')

demands <- map_dfr(1:24, ~make_f(icio=iciot[[.x]], selection=country) %>% 
                          agrega(., axa=2, what='sectoare', selection=country) %>%
                          colSums() %>% as_tibble(rownames='name') %>%
                          mutate(year=.x+1994)) %>%
           rename(demand=value)
message(country)
  x <- rez %>% inner_join(demands)
return(x)
}


sda <- function(icio1, icio2, country, what="M"){
  L1 <- make_L(make_A(icio=icio1, selection=country))
  L2 <- make_L(make_A(icio=icio2, selection=country))
  
  if(what=="M"){
    M1 <- make_M(icio=icio1, selection=country)
    M2 <- make_M(icio=icio2, selection=country)
    }
  else{ 
    M1 <- make_V(icio=icio1, selection=country)
    M2 <- make_V(icio=icio2, selection=country)
    }
  f1 <- make_f(icio=icio1, selection=country)
  f2 <- make_f(icio=icio2, selection=country)
  
  dx1 <- 0.5*(M2-M1)%*%(L1%*%f1+L2%*%f2)
  dx2 <- 0.5*((M1%*%(L2-L1)%*%f2)+(M2%*%(L2-L1)%*%f1))
  dx3 <- 0.5*(M1%*%L1 + M2%*%L2)%*%(f2-f1)
  return(list("dM"=dx1, "dL"=dx2, "df"=dx3))
}


sda_tara <- function(iciot, tara, what="M"){
  imp <- list(23)
  for(i in 1:23){
    imp[[i]] <- sda(icio1=iciot[[i]], icio2=iciot[[(i+1)]], what=what, country = tara) %>% 
      map(., ~agrega2(mat=.x, by='sectoare', selection=tara))
  }
  message(tara)
  return(imp)
}

citeste <- function(year){
  here::here("data", list.files("data", paste0("ICIO2021_", year, ".csv")))
}

#analize==========================================================

contribuţii <- function(country){
bla <- importuri_tara(country = country) 

bla %>% filter(year!=2008) %>%
  group_by(year) %>% mutate(gdp=sum(demand-imp)) %>% ungroup() %>% 
  mutate(adj = demand - imp) %>% mutate(gdp_d=adj/gdp) %>%
  group_by(name) %>% mutate(ratio = 100*((lead(adj)-adj)/adj)*gdp_d) %>%
  drop_na() %>% filter(!is.infinite(ratio)) %>% 
  ungroup() %>% 
  mutate(year = if_else(year<2008, "1995-2007","2009-2018")) %>% 
  group_by(year, name) %>% summarise(across(imp:ratio, ~mean(.,na.rm=T))) %>% 
  ungroup() %>% 
  select(year, name, ratio) %>% spread(year, ratio) %>% 
  mutate(country=country)
}

#=========================regimuri de creştere============================
df <- unique(gsub("_.+$", "", icio$prod))[1:66] %>%
      future_map_dfr(., contribuţii)

nivele <- c("Continental Europe", "LMEs", "Nordic", "Mediterranean/mixed", 
            "Central and Eastern European", "Baltics", "Advanced Asian", 
            "Other Western", "South-East Asian", "China", "India", "African", 
            "Latin and Central American", "Natural resources-rich")

culori <- rep(c("black", "gray50"), 7) 
names(culori) <- nivele


png("contributie.png", width=3000, height = 3200, res=350)
df %>% #group_by(country) %>% 
  #mutate(gdp0 = sum(`1995-2007`), gdp1 = sum(`2009-2018`)) %>%
  #mutate(across(`1995-2007`:`2009-2018`, ~100*./sum(.))) %>% ungroup() %>%
  #filter(gdp1>0|gdp0>0) %>%
  mutate(country=countrycode::countrycode(sourcevar = df$country, 
                                               origin='iso3c', 
                                               destination = 'country.name')) %>%
  mutate_at(2:3, ~if_else(.<0, 0, .)) %>% 
  mutate_at(2:3, ~if_else(.>100, 100, .)) %>%
  inner_join(read_csv(here::here("data/clasificare.csv"))) %>% drop_na() %>%
  group_by(country) %>% 
  mutate(what1=name[which.max(`1995-2007`)], 
         what2 = name[which.max(`2009-2018`)]) %>% 
  filter(name%in%"Exporturi") %>% ungroup() %>% 
  select(-contains("what")) %>% mutate(delta=`2009-2018`-`1995-2007`) %>%
  mutate(country = as.factor(country), clasificare=as.factor(clasificare)) %>% 
  mutate(clasificare = ordered(clasificare, levels=nivele)) %>%
  arrange(clasificare, country) %>% mutate(index=row_number()) %>%
  mutate(culoare = unname(culori[clasificare])) %>% 
  mutate(country = glue::glue("<i style='color:{culoare}'>{country}</i>")) %>%
  mutate(country = fct_reorder(country, index, .desc = TRUE)) %>%
  ggplot(aes(x=country, y=`1995-2007`)) + 
  geom_hline(yintercept = seq(31, 40, length.out=100), colour='grey98') +
  geom_hline(yintercept=seq(40, 54, length.out=100), colour='grey95') +
  geom_hline(yintercept=c(31, 40, 54), linetype=2) +
  geom_point(colour='#EEA100', size=2) + 
  geom_point(aes(x=country, y=`2009-2018`), size=2, colour='royalblue3') + 
  geom_segment(aes(xend=country, x=country, y=`1995-2007`, yend=`2009-2018`,
                   colour=culoare), show.legend = FALSE,
               arrow=arrow(length=unit(0.25, "cm"))) + coord_flip() +
  scale_y_continuous(labels=scales::percent_format(scale=1)) +
  labs(title="Import-adjusted export contribution to growth", 
       subtitle="Source: OECD-ICIO, 2021", x=NULL, y=NULL) +
  scale_colour_manual(values=culori[1:2]) +
  scale_size(range=c(0.1, 5)) +
  facet_grid(rows=vars(clasificare), scales='free_y', space='free_y', switch = 'y') +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size=7), strip.placement = 'outside',
        strip.text.y.left = element_text(angle=0, size=8),
        panel.spacing = unit(0, 'cm'),
        axis.text.y.left = ggtext::element_markdown(),
        strip.text= element_text(family = "Roboto", face = 2),
        plot.title.position = 'plot') 
dev.off()


importuri_tara_sectoare <- function(icio, selection){
va <- make_V(icio=icio, selection=selection)%*%make_L(make_A(icio=icio, selection=selection))%*%make_f(icio=icio, selection=selection)
va <- agrega(va, axa=2, what='sectoare', selection=selection) %>% colSums()

df <- va %>% enframe() %>% rename(prod=name, cerere=value) %>%
      mutate(prod=paste0("D", prod)) %>% 
      left_join(read_csv(here::here("data", 'coduri.csv')))
return(df)
}

#============================================================
#STRUCTURAL DECOMPOSITION ANALYSIS

sda_contrib <- function(tara, what='V', agregare='categorii'){
romania <- sda_tara(iciot, tara=tara, what='V')
x <- map(romania, ~(.x$dL+.x$dM+.x$df) %>% 
      agrega2(mat=., by=agregare, selection=tara) %>% colSums()) %>% 
  tibble(ceva=.) %>% unnest_wider(ceva) %>% 
  mutate(year=1995:2017) %>% relocate(year) %>% 
  gather(2:ncol(.), key='demand', value='valoare') %>%
  mutate(demand = gsub("^.+_", "Internal", demand)) %>% 
  spread(demand, valoare)

if(agregare=="sectoare"){
  return(x %>% mutate(country=tara))
} else{
x<- x %>%
  mutate(Exports = Exporturi + InternalDPABR,
              State = InternalGGFC, 
              Investment = InternalGFCF+InternalINVNT,
              Households = InternalHFCE + InternalNPISH,
              total = Exports + State + Investment + Households) %>%
  select(year, total, Exports, State, Investment, Households) %>%
  gather(3:6, key='component', value='VA') %>% mutate(country=tara)
}
return(x)
}

vax <- future_map(unique(gsub("_.+$", "", icio$prod))[1:66], ~
       sda_contrib(tara=.x, what='V', agregare='categorii')) %>% bind_rows()

vaxd <- future_map(unique(gsub("_.+$", "", icio$prod))[1:66], ~
        sda_contrib(tara=.x, what='V', agregare='sectoare')) %>% bind_rows()


vax %>% 
  mutate(country=countrycode::countrycode(sourcevar=vax$country, 
                                          origin='iso3c', 
                                          destination='country.name')) %>% 
  inner_join(read_csv(here::here("data", 'clasificare.csv'))) %>%
  filter(year!=2008) %>%
  group_by(clasificare, country, component, period=if_else(year<2008, "1995-2007", "2009-2018")) %>%
  mutate(VAR=100*mean(VA, na.rm=TRUE)/abs(mean(total, na.rm=TRUE)), 
         VAp=100*VA/total) %>% ungroup() %>%
  group_by(clasificare, component, period) %>% 
  mutate(VAR_region = 100*mean(VA, na.rm=TRUE)/mean(total, na.rm=TRUE)) %>% 
  ungroup() %>% filter(component=="Exports") %>% 
  mutate_at(vars(VAR:VAp), ~if_else(.<0, 0, .)) %>% 
  mutate_at(vars(VAR:VAp), ~if_else(.>100, 100, .)) %>% drop_na() %>%
  mutate(country = as.factor(country), clasificare=as.factor(clasificare)) %>% 
  mutate(clasificare = ordered(clasificare, levels=nivele)) %>%
  arrange(clasificare, country) %>% group_by(year) %>% 
  mutate(index=row_number()) %>% ungroup() %>%
  mutate(culoare = unname(culori[clasificare])) %>% 
  mutate(country = glue::glue("<i style='color:{culoare}'>{country}</i>")) %>%
  mutate(country = fct_reorder(country, index, .desc = TRUE)) %>%
  select(-year, -total,-VA, -VAp) %>% distinct() %>% 
  pivot_wider(names_from = period, values_from=c(VAR, VAR_region)) %>%
  select(-component) %>% 
  ggplot(aes(x=country, y=`VAR_1995-2007`)) + 
  geom_hline(yintercept = seq(27, 35, length.out=100), colour='grey98') +
  geom_hline(yintercept=seq(35, 46, length.out=100), colour='grey95') +
  geom_hline(yintercept=c(27, 35, 46), linetype=2) +
  geom_point(colour='#EEA100', size=2) + 
  geom_point(aes(x=country, y=`VAR_2009-2018`), size=2, colour='royalblue3') + 
  geom_segment(aes(xend=country, x=country, y=`VAR_1995-2007`, yend=`VAR_2009-2018`,
                   colour=culoare), show.legend = FALSE,
               arrow=arrow(length=unit(0.25, "cm"))) + coord_flip() +
  scale_y_continuous(labels=scales::percent_format(scale=1)) +
  labs(title="Final demand decomposition value-added contributions", 
       subtitle="Source: OECD-ICIO, 2021", x=NULL, y=NULL) +
  scale_colour_manual(values=culori[1:2]) +
  scale_size(range=c(0.1, 5)) +
  facet_grid(rows=vars(clasificare), scales='free_y', space='free_y', switch = 'y') +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size=7), strip.placement = 'outside',
        strip.text.y.left = element_text(angle=0, size=8),
        panel.spacing = unit(0, 'cm'),
        axis.text.y.left = ggtext::element_markdown(),
        strip.text= element_text(family = "Roboto", face = 2),
        plot.title.position = 'plot') 



vax %>% 
  mutate(country=countrycode::countrycode(sourcevar=vax$country, 
                                          origin='iso3c', 
                                          destination='country.name')) %>% 
  inner_join(read_csv(here::here("data", 'clasificare.csv'))) %>%
  filter(year!=2008) %>% 
  group_by(clasificare, year) %>%
  mutate(total=sum(total)) %>%
  ungroup() %>% group_by(clasificare, year, component) %>% 
  mutate(VA=sum(VA)) %>% ungroup() %>%
  group_by(clasificare, component, period=if_else(year<2008, "1996-2007", "2009-2018")) %>%
  summarise(VAR=100*mean(VA, na.rm=TRUE)/abs(mean(total, na.rm=TRUE))) %>% ungroup() %>%
  filter(component=="Exports") %>% drop_na() %>%
  mutate(VAR=if_else(VAR<0, 0, VAR)) %>% 
  mutate(VAR=if_else(VAR>100, 100, VAR)) %>% drop_na() %>%
  mutate(clasificare = ordered(clasificare, levels=nivele)) %>% 
  select(period, VAR, clasificare) %>% distinct() %>% 
  ggplot(aes(x=clasificare, y=VAR, fill=perioad)) + 
  geom_col(aes(fill=period), position='dodge') + coord_flip() +
  scale_fill_tableau() +
  scale_y_log10()


vaxd %>% relocate(country) %>% filter(year!=2008) %>%
  group_by(country, year=if_else(year<2008, "1995-2007", "2009-2018")) %>% 
  summarise(across(2:57, sum)) %>% ungroup() %>% 
  select(1:2, matches("[0-9]+")) %>% 
  gather(3:46, key='component', value='valoare') %>% 
  group_by(country, year) %>% mutate(total=sum(valoare)) %>% 
  filter(grepl("ROU|BGR|HUN|POL|SVK|CZE", country)) %>% 
  mutate(valoare=100*valoare/total) %>% select(-total) %>% 
  ungroup() %>% spread(year, valoare) %>% 
  mutate(delta=`2009-2018`-`1995-2007`) %>% group_by(country) %>% 
  slice_max(abs(delta), n=10) %>% ungroup() %>% 
  rename(prod=component) %>% mutate(prod=paste0("D", prod)) %>% 
  inner_join(read_csv(here::here('data', 'coduri.csv'))) %>% 
  arrange(desc((delta))) %>% select(-prod) %>% 
  relocate(country, sector) %>% group_by(country) %>% 
  gt::gt() %>% gtExtras::gt_theme_538() %>% 
  gt::fmt_percent(2:5, scale_values = FALSE) %>% 
  gt::tab_style(gt::cell_fill(color='#EEA100'), 
                gt::cells_column_labels()) %>% 
  gt::tab_style(gt::cell_fill(color="grey95"), gt::cells_group()) %>% 
  gt::opt_table_font(font='Roboto') %>%
  gt::tab_style(gt::cell_text(weight = 'bold'),gt::cells_column_labels()) %>%
  gt::tab_header(title="Top 10 export sectors by value added generated",
  subtitle="Sorted by the largest absolute difference between periods") %>%
  gt::tab_style(gt::cell_text(weight = 'bold'), gt::cells_title(groups='title')) %>%
  gt::tab_style(gt::cell_text(weight='bold'), gt::cells_group()) %>%
  gt::tab_footnote(footnote="Relative contribution to total VA generated by exports",
                   locations = gt::cells_title(groups='subtitle'))

#==================analize===================================
rez %>% filter(country%in%c("ROU")) %>% 
  filter(year==2018) %>% select(-year) %>%
  gather(2:5, key='FD', value='valoare') %>% 
  mutate(sector=substr(sector, 1, 26)) %>%
  group_by(country) %>%
  mutate(sector = fct_lump(sector, w=valoare, n=30)) %>%
  ungroup() %>%
  mutate(sector = tidytext::reorder_within(sector, valoare, within=country)) %>% 
  ggplot(aes(x=sector, y=valoare, fill="FD")) + 
  geom_col(aes(fill=FD), position='stack') + 
  tidytext::scale_x_reordered() + coord_flip() + 
  facet_wrap(~country, ncol = 3, scales='free') + 
  labs(title = "Driverii importurilor după categorii de cerere finală", 
       x="Sector", y="milioane USD") +
  ggthemes::scale_fill_tableau()

  
#============================================================
i <- importuri(iciot[[24]], selection='ROU')

agrega2 <- function(mat, by='sectoare', selection){
  tmat <- t(mat)
if(by=='tari'){
  tari <- str_extract(colnames(mat), "^[A-Z0-9]{3}")
  x <- t(aggregate.Matrix(tmat, groupings = tari, fun='sum')) %>% as.matrix()
} else if(by=='sectoare'){
  sectoare <- if_else(grepl(selection, colnames(mat)), colnames(mat), 
                      gsub("^[A-Z0-9]{3}_", "", colnames(mat))) 
  x <- t(aggregate.Matrix(tmat, groupings = sectoare, fun='sum')) %>% as.matrix()
} else {
  categorii <- if_else(grepl(selection, colnames(mat)), colnames(mat), "Exporturi")
  x <- t(aggregate.Matrix(tmat, groupings=categorii, fun='sum')) %>% as.matrix()
}
  return(x)
}