
import Control.Monad (replicateM)
import System.Random.MWC (create, GenIO)
import System.Random.MWC.Distributions (beta, bernoulli, dirichlet)
import Numeric.LinearAlgebra.HMatrix
       (Matrix, matrix, Vector, vector, konst, fromRows, (!), toList, fromList, size)
import qualified Data.Vector.Storable as V

import Utils (multinomial)
import VB (lQ, calculateEQZ, calculateQZ, calculateEQmu, calculateQmu, expectationBeta, em)


generateSampleMatrix :: Double  -- ^ alpha parameter for Beta distribution
                     -> Double  -- ^ beta parameter for Beta distribution
                     -> Vector Double  -- ^ alpha parameter for Dirichlet distr
                     -> Int  -- ^ number of samples
                     -> Int  -- ^ vector size
                     -> Int  -- ^ Number of mixins
                     -> IO (Matrix Double, Vector Double)
generateSampleMatrix a b alpha n d k = do
    gen <- create
    pivT <- dirichlet (toList alpha) gen
    let piv = fromList pivT
    muK <-
        replicateM
            k
            (replicateM
                 d
                 (beta a b gen) >>=
             \m ->
                  return (vector m))
    let mu = fromRows muK
    genData <-
        replicateM
            n
            (generateVectors piv mu gen)
    return (fromRows genData, piv)
    where generateVectors :: Vector Double
                          -> Matrix Double
                          -> GenIO
                          -> IO (Vector Double)
          generateVectors piv mu gen = do
              zn <- multinomial piv gen
              let Just j = V.findIndex (== 1) zn
              xn <-
                  mapM
                      (\i ->
                            bernoulli
                                ((mu ! j) !
                                 i)
                                gen >>=
                            \b ->
                                 if b
                                     then return 1.0
                                     else return 0.0)
                      [0 .. d - 1]
              return (vector xn)

-- | L(q) function must increase
testLQ = do
  let k = 6
      n = 1000
      d = 3
      a = 0.6
      b = 1
      alpha = konst (1 / 6) k
  (mX, piV) <- generateSampleMatrix a b alpha n d k
  let test_k = 10
      test_alpha = konst 0.1 test_k
      test_a = matrix d [0.5 | _ <- [1..test_k*d]]
      test_b = matrix d [0.5 | _ <- [1..test_k*d]]
  gen <- create
  pivT <- dirichlet (toList test_alpha) gen
  let piv = fromList pivT
  let qZ = fromRows (replicate n piv)
  let eQz = calculateEQZ qZ
  let qMu = (test_a, test_b)
  let eQmu = calculateEQmu qMu
  print (lQ mX eQz eQmu)
  let qZ' = calculateQZ mX eQmu eQz
      eQz' = calculateEQZ qZ'
      qMu' = calculateQmu mX eQz' qMu
      eQmu' = calculateEQmu qMu'
  print (lQ mX eQz' eQmu')
  let qZ'' = calculateQZ mX eQmu' eQz'
      eQz'' = calculateEQZ qZ''
      qMu'' = calculateQmu mX eQz'' qMu'
      eQmu'' = calculateEQmu qMu''
  print (lQ mX eQz'' eQmu'')

-- | EM first step. Can't imagine.
testEM = do
  let k = 6
      n = 1000
      d = 3
      a = 0.6
      b = 1
      alpha = konst (1 / 6) k
  (mX, piV) <- generateSampleMatrix a b alpha n d k
  let test_k = 10
      test_alpha = konst 0.1 test_k
      test_a = 0.5
      test_b = 0.5
  gen <- create
  pivT <- dirichlet (toList test_alpha) gen
  let piv = fromList pivT
  print "Original pi vector"
  print piV
  print "Calculated pi vector"
  print (em mX piv test_alpha test_a test_b)


main :: IO ()
main = do
  putStrLn "Test suite not yet implemented"
  let k = 3
      n = 10
      d = 3
      a = 1
      b = 1
      alpha = konst 0.01 k
  d <- generateSampleMatrix a b alpha n d k
  print d
  putStrLn "Test L(q). Must increase."
  testLQ
  putStrLn "Test EM work or not :)"
  testEM
  return ()

  
