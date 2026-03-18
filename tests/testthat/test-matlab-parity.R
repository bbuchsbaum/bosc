octave_escape_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  gsub("'", "''", path, fixed = TRUE)
}

find_matlab_source <- function() {
  candidates <- c(
    "behavioral-oscillations",
    file.path("..", "behavioral-oscillations"),
    file.path(testthat::test_path("..", ".."), "behavioral-oscillations")
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    if (dir.exists(candidate)) return(candidate)
  }
  NULL
}

skip_if_no_octave_or_matlab <- function() {
  if (nzchar(Sys.which("octave")) == FALSE) {
    skip("Octave is not available.")
  }
  matlab_src <- find_matlab_source()
  if (is.null(matlab_src)) {
    skip("behavioral-oscillations source not found.")
  }
  matlab_src
}

write_octave_shims <- function(dir) {
  writeLines(
    c(
      "function p = normcdf(x, mu, sigma)",
      "  if nargin < 2 || isempty(mu)",
      "    mu = 0;",
      "  end",
      "  if nargin < 3 || isempty(sigma)",
      "    sigma = 1;",
      "  end",
      "  z = (x - mu) ./ (sigma * sqrt(2));",
      "  p = 0.5 * (1 + erf(z));",
      "end"
    ),
    file.path(dir, "normcdf.m")
  )

  writeLines(
    c(
      "function y = normpdf(x, mu, sigma)",
      "  if nargin < 2 || isempty(mu)",
      "    mu = 0;",
      "  end",
      "  if nargin < 3 || isempty(sigma)",
      "    sigma = 1;",
      "  end",
      "  z = (x - mu) ./ sigma;",
      "  y = exp(-0.5 * z.^2) ./ (sigma * sqrt(2*pi));",
      "end"
    ),
    file.path(dir, "normpdf.m")
  )

  writeLines(
    c(
      "function x = norminv(p, mu, sigma)",
      "  if nargin < 2 || isempty(mu)",
      "    mu = 0;",
      "  end",
      "  if nargin < 3 || isempty(sigma)",
      "    sigma = 1;",
      "  end",
      "  p = min(max(p, 0), 1);",
      "  x = mu + sigma .* sqrt(2) .* erfinv(2 .* p - 1);",
      "end"
    ),
    file.path(dir, "norminv.m")
  )

  writeLines(
    c(
      "function [pks, locs] = findpeaks(data, varargin)",
      "  x = data(:);",
      "  n = length(x);",
      "  locs = [];",
      "  for i = 2:(n-1)",
      "    if x(i) > x(i-1) && x(i) >= x(i+1)",
      "      locs(end+1,1) = i; %#ok<AGROW>",
      "    end",
      "  end",
      "  if isempty(locs)",
      "    pks = [];",
      "    return;",
      "  end",
      "  pks = x(locs);",
      "  npeaks = [];",
      "  sortstr = '';",
      "  i = 1;",
      "  while i <= length(varargin)",
      "    key = lower(varargin{i});",
      "    if i + 1 <= length(varargin)",
      "      val = varargin{i+1};",
      "    else",
      "      val = [];",
      "    end",
      "    if strcmp(key, 'npeaks')",
      "      npeaks = val;",
      "    elseif strcmp(key, 'sortstr')",
      "      sortstr = lower(val);",
      "    end",
      "    i = i + 2;",
      "  end",
      "  if strcmp(sortstr, 'descend')",
      "    [pks, ord] = sort(pks, 'descend');",
      "    locs = locs(ord);",
      "  end",
      "  if ~isempty(npeaks)",
      "    k = min(length(pks), npeaks);",
      "    pks = pks(1:k);",
      "    locs = locs(1:k);",
      "  end",
      "end"
    ),
    file.path(dir, "findpeaks.m")
  )
}

run_octave <- function(code_lines, matlab_src) {
  shim_dir <- tempfile("octave-shims-")
  dir.create(shim_dir, recursive = TRUE)
  on.exit(unlink(shim_dir, recursive = TRUE), add = TRUE)
  write_octave_shims(shim_dir)

  script_file <- tempfile(fileext = ".m")
  on.exit(unlink(script_file), add = TRUE)

  script <- c(
    sprintf("addpath('%s');", octave_escape_path(shim_dir)),
    sprintf("addpath('%s');", octave_escape_path(file.path(matlab_src, "func_external"))),
    sprintf("addpath('%s');", octave_escape_path(file.path(matlab_src, "func_internal"))),
    code_lines
  )
  writeLines(script, script_file)

  output <- suppressWarnings(
    system2(
      "octave",
      c("--quiet", "--no-gui", "--no-window-system", script_file),
      stdout = TRUE,
      stderr = TRUE
    )
  )

  output <- output[!grepl("^error: ignoring const execution_exception& while preparing to exit$", output)]
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(paste(output, collapse = "\n"), call. = FALSE)
  }
  output
}

parse_octave_kv <- function(lines) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & grepl("=", lines, fixed = TRUE)]
  out <- list()
  for (line in lines) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    out[[key]] <- value
  }
  out
}

parse_num_vec <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- gsub(";", " ", x, fixed = TRUE)
  tokens <- strsplit(trimws(x), "[,[:space:]]+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  as.numeric(tokens)
}

test_that("circular statistics agree with MATLAB reference", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "alpha = [0; pi/4; pi/2; -pi/3; 0.2];",
      "w = [1; 2; 1; 3; 2];",
      "mu = circ_mean(alpha, w);",
      "r = circ_r(alpha, w, 0.1);",
      "[pval, v] = circ_vtest(alpha, 0.25, w, 0.1);",
      "fprintf('mu=%.17g\\n', mu);",
      "fprintf('r=%.17g\\n', r);",
      "fprintf('pval=%.17g\\n', pval);",
      "fprintf('v=%.17g\\n', v);"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  alpha <- c(0, pi / 4, pi / 2, -pi / 3, 0.2)
  w <- c(1, 2, 1, 3, 2)
  rv <- circ_vtest(alpha, dir = 0.25, w = w, d = 0.1)

  expect_equal(circ_mean(alpha, w = w), as.numeric(kv$mu), tolerance = 1e-12)
  expect_equal(circ_r(alpha, w = w, d = 0.1), as.numeric(kv$r), tolerance = 1e-12)
  expect_equal(rv$pval, as.numeric(kv$pval), tolerance = 1e-12)
  expect_equal(rv$v, as.numeric(kv$v), tolerance = 1e-12)
})

test_that("pairwise phase consistency agrees with MATLAB reference", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "data = [0, 0.1, NaN; pi/4, 0.2, 0.3; -pi/3, 0.2, 0.4; 0.5, 0.1, 0.5];",
      "ppc = pairwisePhaseConsistency1D(data, 1);",
      "fprintf('ppc=%s\\n', mat2str(ppc, 17));"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  data <- matrix(
    c(
      0, 0.1, NA,
      pi / 4, 0.2, 0.3,
      -pi / 3, 0.2, 0.4,
      0.5, 0.1, 0.5
    ),
    nrow = 4,
    byrow = TRUE
  )
  r_ppc <- as.numeric(pairwise_phase_consistency(data, dim = 1))
  m_ppc <- parse_num_vec(kv$ppc)

  expect_equal(r_ppc, m_ppc, tolerance = 1e-12)
})

test_that("U-score agrees with MATLAB reference", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "dat = [0.8, 1.2, 1.5; 0.2, 0.9, -0.1];",
      "ref = zeros(4, 2, 3);",
      "ref(:,1,1) = [0.1; 0.2; 0.3; 0.4];",
      "ref(:,2,1) = [0.0; 0.1; 0.2; 0.3];",
      "ref(:,1,2) = [0.7; 1.0; 1.3; 1.4];",
      "ref(:,2,2) = [0.2; 0.4; 0.6; 0.8];",
      "ref(:,1,3) = [0.9; 1.1; 1.2; 1.4];",
      "ref(:,2,3) = [-0.3; -0.2; -0.1; 0.0];",
      "[Z, sgnf] = UScoreMat(dat, ref, 0.05);",
      "fprintf('Z=%s\\n', mat2str(Z(:), 17));",
      "fprintf('sgnf=%s\\n', mat2str(sgnf(:), 17));"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  dat <- rbind(c(0.8, 1.2, 1.5), c(0.2, 0.9, -0.1))
  ref <- array(0, dim = c(4, 2, 3))
  ref[, 1, 1] <- c(0.1, 0.2, 0.3, 0.4)
  ref[, 2, 1] <- c(0.0, 0.1, 0.2, 0.3)
  ref[, 1, 2] <- c(0.7, 1.0, 1.3, 1.4)
  ref[, 2, 2] <- c(0.2, 0.4, 0.6, 0.8)
  ref[, 1, 3] <- c(0.9, 1.1, 1.2, 1.4)
  ref[, 2, 3] <- c(-0.3, -0.2, -0.1, 0.0)

  u <- u_score_matrix(dat, ref, alpha = 0.05)
  expect_equal(as.numeric(u$Z), parse_num_vec(kv$Z), tolerance = 1e-12)
  expect_equal(as.numeric(u$sgnf), parse_num_vec(kv$sgnf), tolerance = 1e-12)
})

test_that("nonparam_pval uses documented strict-tail counting convention", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "data = [1.2; -0.5; 0.3];",
      "refv = [-1.0, -0.2, 0.1, 0.6, 1.5];",
      "[h, p] = nonParamPVal(data, refv, 0.1);",
      "fprintf('h=%s\\n', mat2str(h, 17));",
      "fprintf('p=%s\\n', mat2str(p, 17));"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  data <- c(1.2, -0.5, 0.3)
  refv <- c(-1.0, -0.2, 0.1, 0.6, 1.5)
  np <- nonparam_pval(data, refv, alpha = 0.1)
  # Strict-tail counting excludes ties at the observed value.
  expect_equal(as.numeric(np$p), c(0.2, 0.2, 0.4), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(as.numeric(np$p), parse_num_vec(kv$p), tolerance = 1e-12)))
  expect_false(isTRUE(all.equal(as.numeric(np$h), parse_num_vec(kv$h), tolerance = 1e-12)))
})

test_that("make_continuous_trace agrees with MATLAB reference", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "cfg = struct();",
      "cfg.dt = 0.01;",
      "cfg.sd_smooth = 0.02;",
      "cfg.removeVal = 0.20;",
      "cfg.warnings = 'off';",
      "data = [0.01, 0.03, 0.05, 0.09, 0.11, 0.11, NaN, 0.20];",
      "[signal, tspan] = makeContinuousTrace(cfg, data);",
      "fprintf('signal=%s\\n', mat2str(signal, 17));",
      "fprintf('tspan=%s\\n', mat2str(tspan, 17));"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  res <- make_continuous_trace(
    c(0.01, 0.03, 0.05, 0.09, 0.11, 0.11, NA, 0.20),
    dt = 0.01,
    sd_smooth = 0.02,
    remove_val = 0.20,
    warn = FALSE
  )

  expect_equal(res$signal, parse_num_vec(kv$signal), tolerance = 1e-10)
  expect_equal(res$tspan, parse_num_vec(kv$tspan), tolerance = 1e-12)
})

test_that("spectral_peak agrees with MATLAB reference", {
  matlab_src <- skip_if_no_octave_or_matlab()

  output <- run_octave(
    c(
      "cfg = struct();",
      "cfg.fs = 200;",
      "cfg.flim = [3, 20];",
      "cfg.fcor = false;",
      "cfg.taper = 'none';",
      "t = (0:999) / 200;",
      "data = sin(2*pi*8*t) + 0.2*sin(2*pi*12*t);",
      "[freq, fxx, corfft] = spectralPeak(cfg, data);",
      "fprintf('freq=%.17g\\n', freq);",
      "fprintf('fxx=%s\\n', mat2str(fxx, 17));",
      "fprintf('corfft=%s\\n', mat2str(corfft, 17));"
    ),
    matlab_src
  )
  kv <- parse_octave_kv(output)

  t <- (0:999) / 200
  data <- sin(2 * pi * 8 * t) + 0.2 * sin(2 * pi * 12 * t)
  sp <- spectral_peak(data, fs = 200, flim = c(3, 20), fcor = FALSE, taper = "none")

  expect_equal(sp$freq, as.numeric(kv$freq), tolerance = 1e-12)
  expect_equal(sp$fxx, parse_num_vec(kv$fxx), tolerance = 1e-12)
  expect_equal(sp$spectrum, parse_num_vec(kv$corfft), tolerance = 1e-10)
})
