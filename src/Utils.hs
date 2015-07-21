
module Utils where

import Data.Function (on)
import Data.List (sortBy)
import System.Random.MWC (GenIO,uniform)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Numeric.LinearAlgebra.HMatrix (Vector)
import qualified Data.Vector.Storable as V


multinomial :: Vector Double  -- ^ Probabilities of each of the p different outcomes. These should sum to 1.
            -> GenIO
            -> IO (Vector Int)
multinomial pvals gen
    | V.sum pvals > 1.01 =
        error "Utils.multinomial: sum of probabilities should sum to 1"
    | otherwise = do
        y <- uniform gen
        let (_,sample,_) = prepvals y
        return sample
    where prepvals :: Double -> (Double, Vector Int, Bool)
          prepvals x = foldr
                  (\(j,p) (s,v,z) ->
                        if z
                            then (s, v, z)
                            else if s + p >= x
                                     then ( s
                                          , v V.//
                                            [(j, 1)]
                                          , True)
                                     else (s + p, v, z))
                  (0.0, V.replicate k 0, False)
                  (sortBy
                       (compare `on` snd)
                       (zip
                            [0 ..]
                            (V.toList pvals)))
          k = V.length pvals
