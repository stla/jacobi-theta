module Math.JacobiTheta
  (
    jtheta1,
    jtheta2,
    jtheta3,
    jtheta4,
    jtheta1Dash 
  )
  where
import Data.Complex ( imagPart, magnitude, realPart, Complex(..) )

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

modulo :: Double -> Int -> Double
modulo a p = 
  let p' = fromIntegral p
  in
  if a > 0 
    then a - fromIntegral(p * floor(a/p'))  
    else a - fromIntegral(p * ceiling(a/p'))

dologtheta4 :: Cplx -> Cplx -> Int -> Int -> Cplx
dologtheta4 z tau passes maxiter = 
  dologtheta3 (z + 0.5) tau (passes+1) maxiter

dologtheta3 :: Cplx -> Cplx -> Int -> Int -> Cplx
dologtheta3 z tau passes maxiterloc
  | realPart tau2 > 0.6  = dologtheta4 z (tau2 - 1) (passes + 1) maxiterloc
  | realPart tau2 < -0.6 = dologtheta4 z (tau2 + 1) (passes + 1) maxiterloc
  | magnitude tau2 < 0.98 && imagPart tau2 < 0.98 = 
      i_ * pi * tauprime * z * z 
      + dologtheta3 (z * tauprime) tauprime (passes + 1) maxiterloc
      - log(sqrt tau2 / sqrt i_) 
  | otherwise = argtheta3 z tau2 0 maxiterloc
    where
      rPtau = realPart tau
      rPtau2 = if rPtau > 0
        then modulo (rPtau + 1) 2 - 1
        else modulo (rPtau - 1) 2 + 1
      tau2 = rPtau2 :+ imagPart tau
      tauprime = -1 / tau2

argtheta3 :: Cplx -> Cplx -> Int -> Int -> Cplx
argtheta3 z tau passes maxiterloc
  | passes > maxiterloc = error "Reached maximal iteration."
  | iPz < -iPtau / 2 = argtheta3 (-zuse) tau (passes + 1) maxiterloc
  | iPz >= iPtau / 2 = 
      -2 * pi * quotient * i_ * zmin 
      + argtheta3 zmin tau (passes + 1) maxiterloc
      - i_ * pi * tau * quotient * quotient
  | otherwise = calctheta3 zuse tau
    where
      iPz = imagPart z
      iPtau = imagPart tau
      zuse = modulo (realPart z) 1 :+ iPz
      quotient = fromInt $ floor(iPz / iPtau + 0.5)
      zmin = zuse - tau * quotient
      fromInt :: Int -> Cplx
      fromInt = fromIntegral

calctheta3 :: Cplx -> Cplx -> Cplx
calctheta3 z tau = 
    go 1 1
    where
      qw :: Int -> Cplx
      qw n = exp(inpi * (taun + 2 * z)) + exp(inpi * (taun - 2 * z))
        where
          n' = fromIntegral n 
          inpi = i_ * n' * pi 
          taun = n' * tau     
      go n res
        | isNaN modulus = error "NaN has occured in the summation."
        | isInfinite modulus = error "Infinity reached in the summation."
--        | modulus == 0 = error "Zero has occured in the summation."
        | n >= 3 && areClose res resnew = log res
        | otherwise = go (n + 1) resnew
          where
            modulus = magnitude res
            resnew = res + qw n

-------------------------------------------------------------------------------
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

funM :: Cplx -> Cplx -> Cplx
funM z tau = i_ * pi * (z + tau/4)

ljtheta1 :: Cplx -> Cplx -> Cplx
ljtheta1 z tau = ljtheta2 (z - 0.5) tau

-- | First Jacobi theta function
jtheta1 ::
     Complex Double -- ^ z
  -> Complex Double -- ^ q, the nome
  -> Complex Double
jtheta1 z q = exp(ljtheta1 (z/pi) tau)
  where
    tau = getTauFromQ q

ljtheta2 :: Cplx -> Cplx -> Cplx
ljtheta2 z tau = 
  funM z tau + dologtheta3 (z + 0.5 * tau) tau 0 1000

-- | Second Jacobi theta function
jtheta2 ::
     Complex Double -- ^ z
  -> Complex Double -- ^ q, the nome
  -> Complex Double
jtheta2 z q = exp(ljtheta2 (z/pi) tau)
  where
    tau = getTauFromQ q

-- | Third Jacobi theta function
jtheta3 ::
     Complex Double -- ^ z
  -> Complex Double -- ^ q, the nome
  -> Complex Double
jtheta3 z q = exp(dologtheta3 (z/pi) tau 0 1000)
  where
    tau = getTauFromQ q

-- | Fourth Jacobi theta function
jtheta4 ::
     Complex Double -- ^ z
  -> Complex Double -- ^ q, the nome
  -> Complex Double
jtheta4 z q = exp(dologtheta4 (z/pi) tau 0 1000)
  where
    tau = getTauFromQ q

-- | Derivative of the first Jacobi theta function
jtheta1Dash :: 
     Complex Double -- ^ z
  -> Complex Double -- ^ q, the nome
  -> Complex Double
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
