module Main where
import Approx ( approx )
import Data.Complex ( Complex(..) )
import Math.JacobiTheta
    ( jtheta1, jtheta2, jtheta3, jtheta4, jtheta1Dash, jtheta1Dash0, jthetaAB )
import           Test.Tasty       (defaultMain, testGroup)
import           Test.Tasty.HUnit (assertEqual, testCase)

i_ :: Complex Double
i_ = 0.0 :+ 1.0

q :: Complex Double 
q = exp (-pi)

q' :: Complex Double 
q' = exp (-pi/10)

q'' :: Complex Double 
q'' = exp (i_ * pi * tau)
  where
    tau = 2.0 :+ 2.0

main :: IO ()
main = defaultMain $
  testGroup "Tests"
  [ testCase "a jtheta1 value" $ do
      let expected = 1.1816128551455719 :+ 0.59589712760417439
          obtained = jtheta1 (1 :+ 1) q
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "another jtheta1 value" $ do
      let expected = 1.75929905417707 :+ 0.0
          obtained = jtheta1 2 q'
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "yet another jtheta1 value" $ do
      let expected = 0.539843563932874 :+ 0.26400643871132
          obtained = jtheta1 (1 :+ 1) q''
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "a jtheta2 value" $ do
      let expected = 0.74328632006610539 :+ (-0.904159309718008)
          obtained = jtheta2 (1 :+ 1) q
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "a jtheta3 value" $ do
      let expected = 0.86456184935441778 :+ (-0.28488586703507289)
          obtained = jtheta3 (1 :+ 1) q
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "a jtheta4 value" $ do
      let expected = 1.1351891564632007 :+ 0.28517396444192509
          obtained = jtheta4 (1 :+ 1) q
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "a jtheta1Dash value" $ do
      let expected = 0.81117649363854416 :+ (-0.89452803853474627)
          obtained = jtheta1Dash (1 :+ 1) q
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "Jacobi identity" $ do
      let q''' = 0.556 :+ 0.283
          theta2 = jtheta2 0 q'''
          theta3 = jtheta3 0 q'''
          theta4 = jtheta4 0 q'''
          expected = theta3**4
          obtained = theta2**4 + theta4**4
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "Edge case" $ do
      let tau = 0.7792256 :+ 1.0e-7
          q'''' = exp(i_ * pi * tau)
          obtained = jtheta2 0 q''''
          expected = 27.746815969548447 :+ 31.241216782108797
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "jtheta1Dash at 0" $ do
      let q''' = 0.556 :+ 0.283
          expected = jtheta1Dash0 q'''
          obtained = jtheta1Dash 0 q'''
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected),

    testCase "A formula involving some jthetaAB at z=0" $ do
      let q''' = 0.556 :+ 0.283
          expected = (jthetaAB (1/6) 0.5 0 q)**3
          obtained = (jthetaAB (1/6) (1/6) 0 q)**3 + (jthetaAB (1/6) (5/6) 0 q)**3
      assertEqual ""
        (approx 10 obtained)
        (approx 10 expected)

  ]
