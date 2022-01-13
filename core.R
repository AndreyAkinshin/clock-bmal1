library(deSolve)

simulate <- function(context) {
  with(context, {
    G <- function(A, B, s, x) { A * x^s / (B + x^s) }
    dG <- function(A, B, s, x) { A * B * s * x^(s-1) / (B + x^s)^2 }
    L <- function(a, b, s, x) { a / (b + x^s) }
    dL <- function(a, b, s, x) { - a * s * x^(s - 1) / (b + x^s)^2 }
    
    G1 <- function(u) G(A1, B1, s1, u)
    G2 <- function(x) G(A2, B2, s2, x)
    G3 <- function(x) G(A3, B3, s3, x)
    G4 <- function(x) G(A4, B4, s4, x)
    G5 <- function(b) G(A5, B5, s5, b)
    g1 <- function(w) G(a1, b1, sigma1, w)
    L2 <- function(p) L(a2, b2, sigma2, p)
    L3 <- function(p) L(a3, b3, sigma3, p)
    L4 <- function(p) L(a4, b4, sigma4, p)
    L6 <- function(z) L(a6, b6, sigma6, z)
    
    dG1 <- function(u) dG(A1, B1, s1, u)
    dG2 <- function(x) dG(A2, B2, s2, x)
    dG3 <- function(x) dG(A3, B3, s3, x)
    dG4 <- function(x) dG(A4, B4, s4, x)
    dG5 <- function(b) dG(A5, B5, s5, b)
    dg1 <- function(w) dG(a1, b1, sigma1, w)
    dL2 <- function(p) dL(a2, b2, sigma2, p)
    dL3 <- function(p) dL(a3, b3, sigma3, p)
    dL4 <- function(p) dL(a4, b4, sigma4, p)
    dL6 <- function(z) dL(a6, b6, sigma6, z)
    
    model <- function(t, state, parms) {
      p <- as.numeric(state[1])
      u <- as.numeric(state[2])
      w <- as.numeric(state[3])
      z <- as.numeric(state[4])
      x <- as.numeric(state[5])
      b <- as.numeric(state[6])
      list(c(
        dp = k1 * (G1(u) * g1(w) - p),
        du = k2 * (G2(x) * L2(p) - u),
        dw = k3 * (G3(x) * L3(p) - w),
        dz = k4 * (G4(x) * L4(p) - z),
        dx = k5 * (G5(b) - x),
        db = k6 * (C * L6(z) - b)
      ))
    }
    
    times <- seq(0, T, by = 0.01)
    start <- c(p = p0, u = u0, w = w0, z = z0, x = x0, b = b0)
    traj <- as.data.frame(lsoda(start, times, func = model, parms = 0))

    timelineDf <- rbind(
      data.frame(x = times, y = traj[,2], var = "p", type = "p(t) / PER:CRY"),
      data.frame(x = times, y = traj[,3], var = "u", type = "u(t) / PER"),
      data.frame(x = times, y = traj[,4], var = "w", type = "w(t) / CRY"),
      data.frame(x = times, y = traj[,5], var = "z", type = "z(t) / REV-ERB"),
      data.frame(x = times, y = traj[,6], var = "x", type = "x(t) / CLOCKL:BMAL1"),
      data.frame(x = times, y = traj[,7], var = "b", type = "b(t) / BMAL1")
    )
    timelineDf$var <- factor(timelineDf$var, levels = unique(timelineDf$var))
    timelineDf$type <- factor(timelineDf$type, levels = unique(timelineDf$type))
    
    getStar <- function() {
      tol <- 1e-7
      phi <- function(p) {
        uniroot(function(x) x - G5(C * L6(L4(p) * G4(x))), lower = 0, upper = A5, tol = tol)$root
      }
      psi <- function(x) {
        uniroot(function(p) p - G1(L2(p) * G2(x)) * g1(L3(p) * G3(x)), lower = 0, upper = A1 * a1, tol = tol)$root
      }
      xs <- uniroot(function(x) phi(psi(x)) - x, lower = 0, upper = A1 * a1, tol = tol)$root
      ps <- psi(xs)
      us <- L2(ps) * G2(xs)
      ws <- G3(xs) * L3(ps)
      zs <- G4(xs) * L4(ps)
      bs <- C * L6(zs)
      list(p = ps, u = us, w = ws, z = zs, x = xs, b = bs)
    }
    star <- tryCatch(getStar(), error = function(e) list(p = NA, u = NA, w = NA, z = NA, x = NA, b = NA))
    M <- matrix(c(
      -k1, k1 * dG1(star$u) * g1(star$w), k1 * G1(star$u) * dg1(star$w), 0, 0, 0,
      k2 * dL2(star$p) * G2(star$x), -k2, 0, 0, k2 * L2(star$p) * dG2(star$x), 0,
      k3 * dL3(star$p) * G3(star$x), 0, -k3, 0, k3 * L3(star$p) * dG3(star$x), 0,
      k4 * dL4(star$p) * G4(star$x), 0, 0, -k4, k4 * L4(star$p) * dG4(star$x), 0,
      0, 0, 0, 0, -k5, k5 * dG5(star$b),
      0, 0, 0, C * k6 * dL6(star$z), 0, -k6
    ), nrow = 6, ncol = 6, byrow = TRUE)
    eigenM <- tryCatch(eigen(M), error = function(e) list(values = rep(NA, 6), vectors = matrix(rep(NA, 36), ncol = 6)))
    Ml <- eigenM$values
    Me <- eigenM$vectors
    
    proj <- function(e){
      dotProduct <- rowSums(t(t((traj[,2:7]) - star) * e))
      l <- sqrt(sum(e^2))
      dotProduct / l
    }
    traj$e2 <- proj(Re(Me[,2]))
    traj$e3 <- proj(Im(Me[,2]))
    
    
    list(
      context = context,
      model = model,
      traj = traj,
      timelineDf = timelineDf,
      star = star,
      M = M,
      Ml = Ml,
      Me = Me
    )
  })
}

randomCodeName <- function() {
  consonants <- c("b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "n",
                  "p", "q", "r", "s", "t", "v", "w", "x", "y", "z")
  vowels <- c("a", "e", "i", "o", "u")
  paste0(sample(consonants, 1), sample(vowels, 1),
         sample(consonants, 1), sample(vowels, 1),
         sample(consonants, 1))
}
