module Approx where
import Data.Complex

approx0 :: Int -> Double -> Double
approx0 n x = fromInteger (round $ x * (10^n)) / (10.0^^n)

approx :: Int -> Complex Double -> Complex Double
approx n z = approx0 n (realPart z) :+ approx0 n (imagPart z)
