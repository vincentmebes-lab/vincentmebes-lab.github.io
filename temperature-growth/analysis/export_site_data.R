# Builds temperature-growth/data.js for the website from the replication analysis.
# Run from a folder that holds dell2012_full-5.Rmd and its input data
# (climate_panel.dta, AreaWeightedClimate.dta, the WDI API_NY...csv, cru_country_year.csv):
#   LANG=en_US.UTF-8 Rscript export_site_data.R
# Then copy the resulting data.js next to index.html.
if (!file.exists("ws.RData")) {
  knitr::purl("dell2012_full-5.Rmd", "analysis.R", quiet = TRUE)
  source("analysis.R")
  save.image("ws.RData")
} else {
  load("ws.RData")
}
suppressMessages({library(dplyr); library(fixest); library(jsonlite); library(countrycode)})
r3 <- function(x) round(unname(x), 3)

# Table 2, columns 1-4
mods <- list(m1, m2, m3, m4)
vars <- c("wtem","wtem_poor","wpre","wpre_poor","wtem_hot","wpre_hot")
t2 <- lapply(seq_along(mods), function(i) {
  m <- mods[[i]]; b <- coef(m); s <- se(m); p <- pvalue(m)
  out <- list(col = i, N = nobs(m), r2 = r3(r2(m, "wr2")), coefs = list())
  for (v in vars) if (v %in% names(b)) out$coefs[[v]] <- list(b = r3(b[v]), se = r3(s[v]), p = r3(p[v]))
  if (i > 1) { pe <- poor_effect(m, "wtem"); out$poor <- list(b = r3(pe["effect"]), se = r3(pe["se"]), p = r3(pe["p"])) }
  out
})

# Rich-country standard errors for each sample (the Rmd's run() only keeps the point estimate)
rich_se <- function(df) {
  m <- feols(g ~ tem + tem_poor + pre + pre_poor | iso3 + region^year + poor^year,
             data = df, cluster = ~ parent + region_year)
  c(se = unname(se(m)["tem"]), p = unname(pvalue(m)["tem"]))
}
rs <- rbind(rich_se(A), rich_se(min_years(filter(new, year <= 2003))), rich_se(min_years(new)),
            rich_se(min_years(filter(new, !year %in% 2020:2021))),
            rich_se(min_years(filter(new, year >= 2004), k = 15)))
res$se_rich <- rs[, "se"]; res$p_rich <- rs[, "p"]

# Figure 1 points and per-country yearly series (new data)
nm <- function(x) countrycode(x, "iso3c", "country.name", warn = FALSE)
f1o <- f1 |> left_join(select(info, fips60_06, iso3, poor), by = "fips60_06") |>
  mutate(name = nm(iso3)) |> filter(!is.na(name)) |>
  transmute(name, iso3, poor, gdp = round(lngdp2000, 3), tmin = round(tmin, 2), tmax = round(tmax, 2),
            t50 = round(t50s, 2), t00 = round(t00s, 2))
ser <- new |> filter(!is.na(g), !is.na(tem)) |> group_by(iso3) |>
  mutate(dt = tem - mean(tem)) |> ungroup() |>
  transmute(iso3, year, g = round(g, 2), dt = round(dt, 2))
ctry_list <- info |> mutate(name = nm(iso3)) |> filter(iso3 %in% ser$iso3) |>
  select(iso3, name, poor) |> arrange(name)
series <- lapply(split(ser, ser$iso3), function(d) list(y = d$year, g = d$g, dt = d$dt))
se_o <- data.frame(label = rownames(se_tab), se = r3(se_tab$se), p = r3(se_tab$p), row.names = NULL)

out <- list(
  table2 = t2,
  samples = res |> mutate(across(c(rich, poor, se_poor, p_poor, se_rich, p_rich), r3)),
  rolling = roll |> mutate(across(c(poor, lo, hi), r3)),
  periods = period_effects(m_t) |> mutate(across(c(poor_effect, se), r3)),
  se_options = se_o, lags = lag_tab |> mutate(across(c(cumulative_poor, se, p), r3)),
  iv = list(ols = r3(poor_effect(m_ols_iv_sample, "wtem")), iv = r3(poor_effect(m_iv, "fit_wtem"))),
  desc = desc_grp |> mutate(across(where(is.double), ~ round(.x, 2))),
  resid_sd_temp = r3(desc_var$plus_region_yr[1]), raw_sd_temp = r3(desc_var$raw[1]),
  corr_cru = round(cor(cmp$tem_dm, cmp$atem_dm), 3),
  sens = sens |> mutate(across(c(effect, se), r3)),
  fig1 = f1o, countries = ctry_list, series = series)
writeLines(paste0("window.DJO = ", toJSON(out, auto_unbox = TRUE, digits = NA), ";"), "data.js")
cat("bytes:", file.size("data.js"), " countries:", nrow(ctry_list), "\n")
