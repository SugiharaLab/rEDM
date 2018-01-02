library(rgl)
library(here)

# Define some colors
color_green <- "#00CC00"
color_red <- "#FF0000"
color_blue <- "#0000FF"
color_purple <- "#AA00CC"

# Generate Data
generate_lorenz_attractor <- function(out_file = "model_lorenz.Rdata", 
                                      num_points = 2000, 
                                      down_sample = 1,
                                      dt = 0.01)
{
    # model params
    sigma <- 10
    beta <- 8/3
    rho <- 28
    
    # diff_eq
    f <- function(x)
    {
        return(c(sigma * (x[2] - x[1]), 
                 x[1] * (rho - x[3]) - x[2], 
                 x[1] * x[2] - beta * x[3]))
    }
    
    # initialize data 
    dat <- matrix(20, nrow = num_points, ncol = 3)
    
    # generate data using Runge-Kutta
    for (i in 1:(num_points-1))
    {
        xx <- dat[i,]
        for (k in 1:down_sample)
        {
            kx1 <- f(xx)
            # kx2 <- f(xx + dt/2 * kx1)
            # kx3 <- f(xx + dt/2 * kx2)
            # kx4 <- f(xx + dt * kx3)
            # xx <- xx + (kx1 * dt/6 + kx2 * dt/3 + kx3 * dt/3 + kx4 * dt/6)
            xx <- xx + kx1 * dt
        }
        dat[i + 1,] <- xx
    }
    return(dat)
}

# Rescale vector and add buffer space to boundaries
rescale <- function(x, buffer = 0.15, ...)
{
    return((x - min(x, ...)) / (max(x, ...) - min(x, ...)) * (1 - 2 * buffer) + buffer)
}

# Plot lagged reconstruction of time series `v`
draw_lagged_attractor <- function(v, FUN = lines3d, t = 1220, 
                                  idx = seq(from = 1, to = t), tau = 8, 
                                  spacing = c(0, 0, 0),  ...)
{
    FUN <- match.fun(FUN)
    xx <- v[idx + 2 * tau]
    yy <- v[idx + tau]
    zz <- v[idx]
    
    FUN(xx + spacing[1], 
        yy + spacing[2], 
        zz + spacing[3], ...)
    
    return(c(tail(xx, 1) + spacing[1], 
             tail(yy, 1) + spacing[2], 
             tail(zz, 1) + spacing[3]))
}

# Plot Lorenz attractor with projected time series
make_figure_1 <- function(dat, t = 1120, future_len = 40, 
                          rescale_buffer = 0.05, 
                          userMatrix = matrix(c(0.0001, 0.2484, 0.9687, 0.05, 
                                                0.9417, -0.3258, 0.0835, 0, 
                                                0.3367, 0.9122, -0.2339, 0, 
                                                0, 0, 0, 1), nrow = 4, byrow = TRUE),
                          FOV = 15, zoom = 0.7, cex = 1.75, scale = c(1, 1, 1), 
                          windowRect = c(50, 50, 700, 500),
                          plot_file = NULL)
{
    par3d(windowRect = windowRect, scale = scale, 
          userMatrix = userMatrix, 
          zoom = zoom, FOV = FOV, cex = cex)
    
    idx <- seq(from = 1, to = t + future_len)
    x <- rescale(dat[idx, 1], rescale_buffer)
    y <- rescale(-dat[idx, 2], rescale_buffer)
    z <- rescale(dat[idx, 3], rescale_buffer)
    
    past_t <- seq(from = 1, to = t)
    future_t <- seq(from = t, length.out = future_len)
    if (length(future_t) %% 2 == 1)
        future_t <- future_t[-1]
    plot3d(x[past_t], y[past_t], z[past_t], type = "l", lwd = 2, 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    segments3d(x[future_t], y[future_t], z[future_t], lwd = 1)
    
    # draw axes
    lines3d(c(0, 1), 0, 0, lwd = 2)
    text3d(0.75, 0, -0.05, "x", 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(0, c(0, 1), 0, lwd = 2)
    text3d(0, 0.5, -0.1, "y", 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(0, 1, c(0, 1), lwd = 2)
    text3d(0, 1.15, 0.6, "z", 
           family = "serif", font = 3, usePlotmath = FALSE)
    
    # draw projected time series
    theta_x <- -pi/2 - 0.2
    ts_path <- seq(from = 1.5, to = 0, length.out = t)
    lines3d(x[past_t], 
            ts_path * cos(theta_x), 
            ts_path * sin(theta_x), 
            color = color_red, lwd = 2)
    pos_text <- 1.2
    text3d(-0.05, pos_text * cos(theta_x), pos_text * sin(theta_x), "time")
    
    # draw projection from attractor to axis
    xx <- x[t]
    yy <- y[t]
    zz <- z[t]
    points3d(xx, yy, zz, color = "black", size = 8)
    lines3d(c(xx, xx, xx), c(0, yy, yy), c(0, 0, zz), color = color_red, lwd = 1.5)
    
    par3d(windowRect = windowRect, scale = scale, 
          userMatrix = userMatrix, 
          zoom = zoom, FOV = FOV, cex = cex)
    
    if (!is.null(plot_file))
        rgl.postscript(plot_file, fmt = "pdf")
    
    return()
}

# Plot lags of time series x and associated lagged attractor
make_figure_2 <- function(dat, t = 1130, tau = 8, 
                          rescale_buffer = 0.05,
                          userMatrix = matrix(c(0.7570, 0.6533, 0.0131, 0,
                                                -0.1319, 0.1332, 0.9823, 0,
                                                0.6399, -0.7453, 0.1870, 0,
                                                0, 0, 0, 1), nrow = 4, byrow = TRUE), 
                          FOV = 10, zoom = 1, cex = 1.5, scale = c(1, 1, 1), 
                          windowRect = c(50, 50, 700, 500), 
                          plot_file = NULL)
{
    # set up plot
    par3d(windowRect = windowRect)
    #mfrow3d(1, 2)
    layout3d(matrix(1:2, nrow = 1), widths = c(4, 3))
    next3d()
    
    # plot left-side panel
    x <- rescale(dat[, 1], buffer = rescale_buffer)
    draw_lagged_ts <- function(v, t = t, point_t = t, 
                               spacing = 0, 
                               len_before = 420, 
                               len_after = 280, 
                               col = "#EE0000", 
                               pt_col = "#EE00EE", 
                               line_col = "#000088", 
                               lab = "x")
    {
        idx_before <- seq(to = t, from = t - len_before + 1)
        idx_after <- seq(from = t, to = t + len_after)
        
        # time series
        lines3d((idx_before - (t - len_before)) / (len_before + len_after), 
                v[idx_before] / 4 + spacing, 0, lwd = 2, col = col)

        lines3d((idx_after - (t - len_before)) / (len_before + len_after), 
                   v[idx_after] / 4 + spacing, 0, col = col, lwd = 0.5)
        
        # point and line to axis
        points3d((point_t - (t - len_before)) / (len_before + len_after), v[point_t] / 4 + spacing, 0, 
                 col = pt_col, size = 5)
        lines3d((c(t, point_t) - (t - len_before)) / (len_before + len_after), 
                rep.int(v[point_t] / 4 + spacing, 2), c(0, 0), lwd = 1.5, col = line_col)
        text3d(((t + point_t) / 2 - (t - len_before)) / (len_before + len_after), 
               v[point_t] / 4 * 0.85 + spacing, 0, "tau", cex = 1)
        
        # box and labels
        lines3d(c(0, 0, 1, 1, 0), c(1, 0, 0, 1, 1) / 4 + spacing, c(0, 0, 0, 0, 0))
        text3d(-0.2, 0.5 / 4 + spacing, 0, lab, adj = 0, col = col)
    }
    draw_lagged_ts(x, t + 2 * tau, t + 2 * tau, 0.3, lab = "x_t")
    draw_lagged_ts(x, t - 10 * tau, t + 2 * tau, 0, lab = "x_t-tau", col = "#BB0000")
    draw_lagged_ts(x, t - 22 * tau, t + 2 * tau, -0.3, lab = "x_t-2tau", col = "#880000")
    lines3d(c(0.6, 0.6), c(-0.3, 0.55), c(0, 0), lwd = 2)
    text3d(0.6, -0.34, 0, "t", cex = 1)
    par3d(userMatrix = matrix(c(1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1, 0,
                                0, 0, 0, 1), nrow = 4, byrow = TRUE), 
          FOV = 10, zoom = 0.9)
    
    # plot right-side panel
    next3d()
    par3d(userMatrix = userMatrix, 
          zoom = zoom, FOV = FOV, cex = cex)
    
    x <- rescale(dat[, 1], buffer = rescale_buffer)
    m_x_t <- draw_lagged_attractor(x, idx = seq(from = 1, to = t), lwd = 2)
    points3d(m_x_t[1], m_x_t[2], m_x_t[3], color = "black", size = 8)
    draw_lagged_attractor(x, idx = seq(from = t + 1, length.out = 40), lwd = 0.5)
    lines3d(c(0, 1), 0, 0, lwd = 2, col = "#EE0000")
    text3d(0.5, 0, -0.05, "x_t", col = "#EE0000", 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(1, c(0, 1), 0, lwd = 2, col = "#BB0000")
    text3d(1, 0.75, -0.1, "x_t-tau", col = "#BB0000", 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(0, 0, c(0, 1), lwd = 2, col = "#880000")
    text3d(0.1, 0.1, 1.0, "x_t-2tau", col = "#880000", 
           family = "serif", font = 3, usePlotmath = FALSE)
    text3d(0.2, 0.5, 0.9, "M_x", cex = 1.75)
    
    if (!is.null(plot_file))
        rgl.postscript(plot_file, fmt = "pdf")
    
    return()
}

# Plot reconstructed attractors using lags of x and y
make_figure_3 <- function(dat, t = 1220, 
                          shift = 0.55, 
                          rescale_buffer = 0.05, 
                          userMatrix = matrix(c(0.9809, 0.1941, -0.0105, 0, 
                                                -0.0382, 0.2453, 0.9687, 0, 
                                                0.1906, -0.9498, 0.2481, 0, 
                                                0, 0, 0, 1), nrow = 4, byrow = TRUE),
                          FOV = 10, zoom = 0.7, cex = 1.75, scale = c(1, 1, 1), 
                          windowRect = c(50, 50, 700, 500),
                          plot_file = NULL)
{
    par3d(windowRect = windowRect, scale = scale, 
          userMatrix = userMatrix, 
          zoom = zoom, FOV = FOV, cex = cex)
    
    x <- rescale(dat[, 1], buffer = rescale_buffer)
    y <- rescale(dat[, 2], buffer = rescale_buffer)
    
    plot3d(0, 0, 0, type = "n", 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    
    spacing_x <- c(-shift, 0, 0)
    m_x_t <- draw_lagged_attractor(x, spacing = spacing_x, lwd = 2)
    points3d(m_x_t[1], m_x_t[2], m_x_t[3], color = "black", size = 8)
    draw_lagged_attractor(x, idx = seq(from = t + 1, length.out = 80), 
                          spacing = spacing_x, lwd = 0.5)
    text3d(0.2 + spacing_x[1], 0.5 + spacing_x[2], 1 + spacing_x[3], "M_x", cex = 2)
    
    spacing_y <- c(shift, 0, 0)
    m_y_t <- draw_lagged_attractor(y, spacing = spacing_y, lwd = 2)
    points3d(m_y_t[1], m_y_t[2], m_y_t[3], color = "black", size = 8)
    draw_lagged_attractor(y, idx = seq(from = t + 1, length.out = 80), 
                          spacing = spacing_y, lwd = 0.5)
    text3d(1 + spacing_y[1], 0.5 + spacing_y[2], 1 + spacing_y[3], "M_y", cex = 2)
    
    lines3d(rbind(m_x_t, m_y_t), color = color_blue)
    
    if (!is.null(plot_file))
        rgl.postscript(plot_file, fmt = "pdf")
    
    return()
}

if(TRUE)
{
    dat <- generate_lorenz_attractor()
    

    fig_1_file <- here("vignettes", "vignette_figs", "figure_1.pdf")
    make_figure_1(dat, plot_file = fig_1_file)
    
    open3d()

    fig_2_file <- here("vignettes", "vignette_figs", "figure_2.pdf")
    make_figure_2(dat, plot_file = fig_2_file)
    
    open3d()
    
    fig_3_file <- here("vignettes", "vignette_figs", "figure_3.pdf")
    make_figure_3(dat, plot_file = fig_3_file)
}

# 

