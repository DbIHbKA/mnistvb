
module VB where

import Numeric.SpecFunctions (digamma)
import Numeric.LinearAlgebra.HMatrix
       (Matrix, size, (!), fromRows, scalar, sumElements, toRows, reshape,
        flatten, tr, dot, (<>), ( #> ), flatten, asRow, matrix, fromList, Vector, cmap)
import Numeric.LinearAlgebra.Devel (zipVectorWith)


data LogType
    = Log  -- ^ log x
    | InvLog  -- ^ log (1 - x)


type ExpectationQZ = Matrix Double
type ExpectationQmu = (Matrix Double, Matrix Double)  -- ^ Matrix Log and Matrix InvLog
type QZ = Matrix Double
type Qmu = (Matrix Double, Matrix Double)  -- ^ Matrix A and matrix B

-- | Moments of logarithmically transformed random variables
expectationBeta :: Double  -- ^ alpha
                -> Double -- ^ beta
                -> LogType  -- ^ log x or log (1 - x)
                -> Double
expectationBeta alpha beta (Log) = digamma alpha -
    digamma (alpha + beta)
expectationBeta alpha beta (InvLog) = digamma beta -
    digamma (alpha + beta)


em :: Matrix Double
   -> Vector Double
   -> Vector Double
   -> Double
   -> Double
   -> Vector Double
em = go 0 (-1)
  where
    go lQthetaNew lQthetaOld matrixX vectorPi vectorAlpha a b
        | abs (lQthetaNew - lQthetaOld) <
              0.1 =
            vectorPi
        | otherwise =
            let (eQz,eQmu) = eStep matrixX initQZ initQmu
                vectorPiNew = mStep vectorPi vectorAlpha eQz
            in go
                   (lQTheta matrixX vectorPiNew vectorAlpha eQz eQmu)
                   lQthetaNew
                   matrixX
                   vectorPiNew
                   vectorAlpha
                   a
                   b
        where initQZ = fromRows (replicate nRow vectorPi)
              initQmu = (matrixAB a, matrixAB b)
              matrixAB v = matrix
                      nColumn
                      [v | _ <-
                              [1 .. nMixtures * nColumn]]
              (nRow,nColumn) = size matrixX


eStep :: Matrix Double -> QZ -> Qmu -> (ExpectationQZ, ExpectationQmu)
eStep matrixX qZ qMu = let eQmu = calculateEQmu qMu
                           eQz = calculateEQZ qZ
                           qZ' = calculateQZ matrixX eQmu qZ
                           eQz' = calculateEQZ qZ'
                           qMu' = calculateQmu matrixX eQz' qMu
                           eQmu' = calculateEQmu qMu'
    in go
           (lQ matrixX eQz' eQmu')
           (lQ matrixX eQz eQmu)
           qZ'
           eQz'
           qMu'
           eQmu'
    where go lQnew lQold qZn eQzn qMun eQmun
              | abs (lQnew - lQold) <
                    0.1 =
                  (eQzn, eQmun)
              | otherwise =
                  let qZn' = calculateQZ matrixX eQmun qZn
                      eQzn' = calculateEQZ qZn'
                      qMun' = calculateQmu matrixX eQzn' qMun
                      eQmun' = calculateEQmu qMun'
                  in go (lQ matrixX eQzn' eQmun') lQnew qZn' eQzn' qMun' eQmun'


lQ :: Matrix Double -> ExpectationQZ -> ExpectationQmu -> Double
lQ matrixX eQz (eQmuM,eQmuMi) = sum
        (map
             (\n ->
                   (eQz ! n) `dot`
                   flatten
                       (asRow (matrixX ! n) <>
                        teQmuM +
                        asRow
                             (1 -
                              (matrixX ! n)) <>
                        teQmuMi))
             [0 .. nRow - 1])
    where teQmuM = tr eQmuM
          teQmuMi = tr eQmuMi
          (nRow,_) = size matrixX

calculateEQZ :: QZ -> ExpectationQZ
calculateEQZ qZ = fromRows
        (map
             (\row ->
                   row /
                   scalar (sumElements row))
             (toRows qZ))

calculateQZ :: Matrix Double -> ExpectationQmu -> QZ -> QZ
calculateQZ matrixX (eQmuM,eQmuMi) qZ = qZ + matrixX <> tr eQmuM +
    (1 - matrixX) <>
    tr eQmuMi


calculateEQmu :: Qmu -> ExpectationQmu
calculateEQmu (aM, bM) = (reshape d qMuM, reshape d qMuMi)
  where
    aV = flatten aM
    bV = flatten bM
    qMuM = zipVectorWith (\a b -> expectationBeta a b Log) aV bV
    qMuMi = zipVectorWith (\a b -> expectationBeta a b InvLog) aV bV
    (k, d) = size aM

calculateQmu :: Matrix Double -> ExpectationQZ -> Qmu -> Qmu
calculateQmu matrixX eQz (aOld,bOld) = ( aOld - 1 + tr eQz <> matrixX
                                       , bOld - 1 +
                                         tr eQz <>
                                         (1 - matrixX))


lQTheta :: Matrix Double
        -> Vector Double
        -> Vector Double
        -> ExpectationQZ
        -> ExpectationQmu
        -> Double
lQTheta matrixX vectorPi vectorAlpha eQz eQmu = lQ matrixX eQz eQmu +
    vectorAlpha `dot`
    cmap log vectorPi

mStep :: Vector Double -> Vector Double -> ExpectationQZ -> Vector Double
mStep vectorPi vectorAlpha eQz = let vectorNewPi = fromList
                                             (map
                                                  (\k ->
                                                        -(vectorPi ! k) *
                                                         sumElements
                                                              (eqZt ! k) /
                                                         (vectorAlpha ! k - 1))
                                                  [0 .. nMixtures - 1])
    in norm vectorNewPi
    where eqZt = tr eQz
          norm v = v /
              scalar (sumElements v)



nMixtures :: Int
nMixtures = 10
