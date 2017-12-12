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

make_attractor_plot <- function(x, y, z, 
                                t = floor(length(x) / 2), 
                                rescale_buffer = 0.05, 
                                userMatrix = matrix(c(0.0001, 0.2484, 0.9687, 0.05, 
                                                      0.9417, -0.3258, 0.0835, 0, 
                                                      0.3367, 0.9122, -0.2339, 0, 
                                                      0, 0, 0, 1), nrow = 4, byrow = TRUE),
                                FOV = 15, zoom = 0.7, 
                                scale = c(1, 1, 1), 
                                windowRect = c(50, 50, 1000, 750),
                                plot_file = NULL)
{
    x <- rescale(x, rescale_buffer)
    y <- rescale(-y, rescale_buffer)
    z <- rescale(z, rescale_buffer)
    
    past_t <- seq(from = 1, to = t)
    future_t <- seq(from = t, to = length(x))
    if (length(future_t) %% 2 == 1)
        future_t <- future_t[-1]
    plot3d(x[past_t], y[past_t], z[past_t], type = "l", lwd = 2, 
           xlab = "", ylab = "", zlab = "", axes = FALSE)
    segments3d(x[future_t], y[future_t], z[future_t], lwd = 1)

    # draw axes
    lines3d(c(0, 1), 0, 0, lwd = 2)
    text3d(0.75, 0, -0.05, "x", cex = 1.5, 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(0, c(0, 1), 0, lwd = 2)
    text3d(0, 0.5, -0.1, "y", cex = 1.5, 
           family = "serif", font = 3, usePlotmath = FALSE)
    lines3d(0, 1, c(0, 1), lwd = 2)
    text3d(0, 1.15, 0.6, "z", cex = 1.5, 
           family = "serif", font = 3, usePlotmath = FALSE)
    
    # draw projected time series
    theta_x <- -pi/2 - 0.2
    ts_path <- seq(from = 1.5, to = 0, length.out = t)
    lines3d(x[past_t], 
            ts_path * cos(theta_x), 
            ts_path * sin(theta_x), 
            color = color_red, lwd = 2)
    pos_text <- 1
    text3d(-0.05, pos_text * cos(theta_x), pos_text * sin(theta_x), "time", 
           cex = 1.5)
    arrow_start <- 0.85
    arrow_end <- 0.5
    lines3d(x = c(-0.05, -0.05), 
            y = c(arrow_start * cos(theta_x), arrow_end * cos(theta_x)), 
            z = c(arrow_start * sin(theta_x), arrow_end * sin(theta_x)), 
            lwd = 2)
    # arrow3d(c(-0.05, arrow_start * cos(theta_x), -arrow_start * sin(theta_x)), 
    #         c(-0.05, arrow_end * cos(theta_x), -arrow_end * sin(theta_x)), 
    #         type = "flat", s = 0.15, theta = pi/8)
    
    # draw projection from attractor to axis
    xx <- x[t]
    yy <- y[t]
    zz <- z[t]
    points3d(xx, yy, zz, color = "black", size = 8)
    lines3d(c(xx, xx, xx), c(0, yy, yy), c(0, 0, zz), color = color_red, lwd = 1.5)

    par3d(windowRect = windowRect, scale = scale, 
        userMatrix = userMatrix, 
        zoom = zoom, FOV = FOV)

    if (!is.null(plot_file))
        rgl.postscript(plot_file, fmt = "pdf")
    
    return()
}

fig_1_file <- here("vignettes", "vignette_figs", "figure_1.pdf")
dat <- generate_lorenz_attractor()
make_attractor_plot(dat[1:1160, 1],
                    dat[1:1160, 2],
                    dat[1:1160, 3], 
                    t = 1120, 
                    plot_file = fig_1_file)

# img <- image_read(fig_1_file)
# print(img)
# image_trim(img)
# print(img)
# image_write(img, path = fig_1_file)


