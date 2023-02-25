module Math.JacobiTheta
  (
    jtheta1Dash,
    jtheta1,
    jtheta2,
    jtheta3,
    jtheta4
  )
  where
import Data.Complex

type Cplx = Complex Double

i_ :: Cplx
i_ = 0.0 :+ 1.0

machinePrecision :: Double
machinePrecision = 2**(-52)

areClose :: Cplx -> Cplx -> Bool
areClose z1 z2 = magnitude (z1 - z2) < epsilon * h
  where
    epsilon = 2.0 * machinePrecision
    magn2 = magnitude z2
    h = if magn2 < epsilon then 1.0 else max (magnitude z1) magn2

square :: Cplx -> Cplx
square z = z * z

jtheta1Alt1 :: Cplx -> Cplx -> Cplx
jtheta1Alt1 z q =
  go 0 (0.0 :+ 0.0) 1.0 (1.0 / qsq) 1.0
  where 
    qsq = q * q
    go :: Int -> Cplx -> Cplx -> Cplx -> Cplx -> Cplx
    go n out alt q_2n q_n_np1 
      | n > 3000 = error "Reached 3000 iterations."
      | areClose out outnew = 2.0 * sqrt (sqrt q) * out
      | otherwise = go (n + 1) outnew (-alt) q_2np1 q_np1_np2
        where
          q_2np1 = q_2n * qsq
          q_np1_np2 = q_n_np1 * q_2np1
          n' = fromIntegral n 
          k = 2.0 * n' + 1.0
          outnew = out + alt * q_np1_np2 * sin (k * z) 

-- jtheta1(z, tau) = jtheta1Alt2 (z/pi) (-i_ * tau/pi)
jtheta1Alt2 :: Cplx -> Cplx -> Cplx
jtheta1Alt2 z' t' = 
  let nm = round (0.5 - realPart z') in
  let np = nm + 1 in
  go nm np (0.0 :+ 0.0) (if even np then (-1, 1) else (1, -1)) 
  where
    go :: Int -> Int -> Cplx -> (Cplx, Cplx) -> Cplx
    go nminus nplus series (altm, altp)
      | nplus - nminus > 3000 = error "Reached 3000 iterations."
      | (nplus - nminus > 2) && areClose series newseries = 
          series / sqrt (pi * t')
      | otherwise = go (nminus - 1) (nplus + 1) newseries (-altm, -altp)
        where 
          nminus' = fromIntegral nminus
          nplus' = fromIntegral nplus
          termm = altm * exp (- square (nminus' - 0.5 + z') / t')
          termp = altp * exp (- square (nplus' - 0.5 + z') / t')
          newseries = series + termm + termp

falpha :: Cplx -> Cplx -> Cplx
falpha z tau = 
  sqrt (-i_ * tau) * exp (i_ / tau * z * z / pi)

jtheta1Alt :: Cplx -> Cplx -> Cplx
jtheta1Alt z tau = 
  if imagPart tau > 1.3 
    then
      let w = pi * tau in 
      i_ * jtheta1Alt2 (z / w) (i_ / w) / falpha z tau
    else
      i_ * jtheta1Alt1 (z / tau) (exp (-i_ * pi / tau)) / falpha z tau

tauFromQ :: Cplx -> Cplx
tauFromQ q = -i_ * log q / pi

checkQ :: Cplx -> Cplx
checkQ q
  | magnitude q >= 1 = 
    error "The modulus of the nome must be smaller than one."
  | imagPart q == 0 && realPart q <= 0 = 
    error "If the nome is real, it must be positive."
  | otherwise = q

getTauFromQ :: Cplx -> Cplx
getTauFromQ = tauFromQ . checkQ

expM :: Cplx -> Cplx -> Cplx
expM z tau = exp (i_ * (z + tau * pi/4))

-- | First Jacobi theta function
jtheta1 :: 
     Cplx -- ^ z
  -> Cplx -- ^ q, the nome
  -> Cplx
jtheta1 z q = jtheta1Alt z (getTauFromQ q)

-- | Second Jacobi theta function
jtheta2 :: 
     Cplx -- ^ z
  -> Cplx -- ^ q, the nome
  -> Cplx
jtheta2 z = jtheta1 (z + pi/2)

-- | Third Jacobi theta function
jtheta3 :: 
     Cplx -- ^ z
  -> Cplx -- ^ q, the nome
  -> Cplx
jtheta3 z q = jtheta2 (z - pi/2 * tau) q * expM (-z) tau
  where
    tau = tauFromQ q

-- | Fourth Jacobi theta function
jtheta4 :: 
     Cplx -- ^ z
  -> Cplx -- ^ q, the nome
  -> Cplx
jtheta4 z = jtheta3 (z + pi/2)

-- | Derivative of the first Jacobi theta function
jtheta1Dash :: 
     Cplx -- ^ z
  -> Cplx -- ^ q, the nome
  -> Cplx
jtheta1Dash z q = 
  go 0 (0.0 :+ 0.0) 1.0 (1.0 / qsq) 1.0
  where 
    q' = checkQ q
    qsq = q' * q'
    go :: Int -> Cplx -> Cplx -> Cplx -> Cplx -> Cplx
    go n out alt q_2n q_n_np1 
      | n > 3000 = error "Reached 3000 iterations."
      | areClose out outnew = 2.0 * sqrt (sqrt q) * out
      | otherwise = go (n + 1) outnew (-alt) q_2np1 q_np1_np2
        where
          q_2np1 = q_2n * qsq
          q_np1_np2 = q_n_np1 * q_2np1
          n' = fromIntegral n 
          k = 2.0 * n' + 1.0
          outnew = out + k * alt * q_np1_np2 * cos (k * z) 
