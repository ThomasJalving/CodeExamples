{-# OPTIONS_GHC -Wno-incomplete-patterns #-}
module PingPong.Player.MostBestPlayer (player, collision, plan) where

import Control.Lens
import Data.Geometry
import Data.Ext
import Data.Fixed

import PingPong.Model
import PingPong.Player

import PingPong.Simulation.Collision

import Data.Colour
import Data.Colour.Names
import PingPong.Simulation (collide', modAngle, simpleBallStep, unJust)
import Transformation (reflection)
import Data.Maybe (isJust)
import Data.Tree.Util (root)
import Data.Geometry.Vector (unit)
import Data.Geometry.Vector.VectorFamily(vectorFromListUnsafe)
import Debug.Trace
import Data.Geometry.Matrix
import qualified Data.Geometry.Vector.VectorFamilyPeano as Data.Geometry.Vector
import Data.Geometry.Ball (Ball)
import PingPong.Model (BallState)



player :: Player
player = defaultPlayer
  { name    = "Most Best Player"
  , arm     = bestArm
  , foot    = bestFoot
  , action  = bestAction
  }

gradient :: Float -> Colour Float
gradient x = blend x goldenrod yellow

bestArm :: Arm
bestArm = [ Link  (gradient 0.1) 0.5
          , Joint (gradient 0.2) (-0.5 * Prelude.pi)
          , Link  (gradient 0.3) 0.5
          , Joint (gradient 0.4) (Prelude.pi)
          , Link  (gradient 0.5) 0.5
          , Joint (gradient 0.6) (-Prelude.pi)
          , Link  (gradient 0.7) 0.5
          , Joint (gradient 0.8) (Prelude.pi)
          , Link  (gradient 0.9) 0.5
          , Joint (gradient 0.6) (-0.5 * Prelude.pi)
          , Link  (gradient 0.7) 0.1
          ]

bestFoot :: Float
bestFoot = 1.5

neutralArmAngles :: [Float]
neutralArmAngles = armToAngles bestArm

bestAction :: Float -> (Float, Item) -> BallState -> Arm -> IO Motion
bestAction t (_, Bat Opponent) bs arm = do
  (collisionBallState, tempTime) <- getBallStateAfterCollision t bs
  (highestBallState,targetTime) <- getBallStateHighestPoint tempTime collisionBallState
  armPositionBestAction 20 t targetTime highestBallState arm
bestAction t (_, Table Self) bs arm = do
  (highestBallState,targetTime) <- getBallStateHighestPoint t bs
  armPositionBestAction 20 t targetTime highestBallState arm
bestAction t _ _ arm = return $ returnToNeutralStance t arm

armPositionBestAction :: Float -> Float -> Float -> BallState -> Arm -> IO Motion
armPositionBestAction maxForward currentTime targetTime bs arm = do
  ccdResult <- ccdprecompute (bestFoot, arm) (loc bs, Vector2 (-1) 0, Vector2 0 0)
  let (angles, angulvels) = unJust ccdResult
  let possible = all (\(a, b) -> abs ((a - b)/2) <= targetTime - currentTime + 0.01) (reduceAngleDifferenceToPi (zip (armToAngles arm) angles))
  let (Point2 _ currentY) = loc bs
  if maxForward > 0 && isJust ccdResult && not possible && currentY > 0.5 then armPositionBestAction (maxForward - 1) currentTime (targetTime + 0.02) (simpleBallStep 0.02 bs) arm
    else if isJust ccdResult
        then return $ moveToStance currentTime targetTime arm angles angulvels
        else return $ returnToNeutralStance currentTime arm


simpleBallStepBackwards :: Float -> BallState -> BallState
simpleBallStepBackwards f st =
  let decay = (1 - 0.05) ** (-f) in
    let newdir = decay *^ dir st ^+^ 2 * (-f) *^ down
    in st { loc = loc st .+^ (-f) *^ newdir
          , dir = newdir
          }

down :: Vector 2 Float
down = Vector2 0 (-1)

getBallStateAfterCollision :: Float -> BallState -> IO (BallState,Float)
getBallStateAfterCollision t st@(BallState _ (Vector2 _ y)) | y <= 0  = do

  let nextState = simpleBallStep 0.02 st
      tableSegment = ClosedLineSegment (Point2 0 0.5 :+ ()) (Point2 1 0.5 :+ ())
  maybeCollCheck <- collision (0, loc st, tableSegment) (0.02, loc nextState , tableSegment)
  let (t, pos, dir) = unJust maybeCollCheck
      finalState = if isJust maybeCollCheck then simpleBallStep (0.02 - t) (BallState pos dir) else nextState
  getBallStateAfterCollision (t + 0.02) finalState
                                                          | otherwise  = return (st,t)

getBallStateHighestPoint :: Float -> BallState -> IO (BallState, Float)
getBallStateHighestPoint t st@(BallState _ (Vector2 _ y)) | y >= 0  = getBallStateHighestPoint (t + 0.02) (simpleBallStep 0.02 st)
                                                        | otherwise  = return (st , t)

moveToStance :: Float -> Float -> Arm -> [Angle] -> [AngulVel] -> Motion
moveToStance startTime targetTime arm targetAngles targetVelocities = let
  tempTargetAngles = convertAngles targetAngles
  in let armAngles = armToAngles arm
  in let (realArmAngles, realTargetAngles) = unzip $ reduceAngleDifferenceToPi (zip armAngles tempTargetAngles)
  in if targetTime - startTime <= 0.1 
    then zipWith3 (calculateAngleVelocity3 startTime targetTime) realArmAngles realTargetAngles targetVelocities
    else  map (\a -> a/abs a * 2) $ zipWith (-) realTargetAngles realArmAngles

calculateAngleVelocity3 :: Float -> Float -> Angle -> Angle -> AngulVel -> AngulVel
calculateAngleVelocity3 startTime targetTime startAngle targetAngle targetVelocity = let
  startMatrix = [[startTime*startTime, startTime, 1],
                 [targetTime*targetTime, targetTime, 1],
                 [2 * targetTime, 1, 0]]
  in let resultMatrix = [[startAngle], [targetAngle], [targetVelocity]]
  in let [[a],[b],[c]] = multiplyMatrix (inverseMatrix33 startMatrix) resultMatrix
  in (2*a*startTime) + b

calculateAngleVelocity4 :: Float -> Float -> Angle -> Angle -> AngulVel -> AngulVel
calculateAngleVelocity4 startTime targetTime startAngle targetAngle targetVelocity = let
  startMatrix = [[startTime*startTime*startTime, startTime*startTime, startTime, 1],
                 [targetTime*targetTime*targetTime, targetTime*targetTime, targetTime, 1],
                 [3 * 4, 2 * 2, 1, 0],
                 [3 * (targetTime * targetTime), 2 * targetTime, 1, 0]]
  in let resultMatrix = [[startAngle], [targetAngle], [2], [targetVelocity]]
  in let [[a],[b],[c],[d]] = multiplyMatrix (inverseMatrix44 startMatrix) resultMatrix
  in (3*a*(startTime * startTime) + (2 * b * startTime) + startTime)

convertAngles :: [Angle] -> [Angle]
convertAngles = map (\angle -> if angle < 0 then angle - ((2 * Prelude.pi) * fromInteger(fst (divMod' angle (2*Prelude.pi)))) else angle)

reduceAngleDifferenceToPi :: [(Angle,Angle)] -> [(Angle,Angle)]
reduceAngleDifferenceToPi = map (\tuple -> if abs (fst tuple - snd tuple) > Prelude.pi 
  then (if fst tuple < snd tuple then (fst tuple + (2*Prelude.pi), snd tuple) else (fst tuple, snd tuple + (2*Prelude.pi))) 
  else tuple)

returnToNeutralStance :: Float -> Arm -> Motion
returnToNeutralStance startTime arm = moveToStance startTime (startTime + 0.2) arm neutralArmAngles [0,0,0,0,0]

inverseMatrix33 :: CustomMatrix -> CustomMatrix
inverseMatrix33 [[a, b, c],
                [d, e, f],
                [g, h, i]] = let matrix = Matrix $ Vector3 (Vector3 a b c) (Vector3 d e f) (Vector3 g h i)
                in let inverseMatrix = inverse' matrix
                in let Matrix (Vector3 (Vector3 a' b' c') (Vector3 d' e' f') (Vector3 g' h' i')) = inverseMatrix
                in [[a' ,b' ,c'],
                [d', e', f'],
                [g', h', i']]

inverseMatrix44 :: CustomMatrix -> CustomMatrix
inverseMatrix44 [[a, b, c, d],
                [e, f, g, h],
                [i, j, k, l],
                [m, n, o, p]] = let matrix = Matrix $ Vector4 (Vector4 a b c d) (Vector4 e f g h) (Vector4 i j k l) (Vector4 m n o p)
                in let inverseMatrix = inverse' matrix
                in let Matrix (Vector4 (Vector4 a' b' c' d') (Vector4 e' f' g' h') (Vector4 i' j' k' l') (Vector4 m' n' o' p')) = inverseMatrix
                in [[a', b', c', d'],
                [e', f', g', h'],
                [i', j', k', l'],
                [m', n', o', p']]

collision :: CollisionChecker
collision (time1, point1, segment1) (time2, point2, segment2) = do
  let Point2 xp1 yp1 = point1
      Point2 xp2 yp2 = point2
      Point2 xq1 yq1 = segment1 ^. start . core
      Point2 xq2 yq2 = segment2 ^. start . core
      Point2 xr1 yr1 = segment1 ^. end . core
      Point2 xr2 yr2 = segment2 ^. end . core
      a              = yp1 - yq1
      b              = yp2 - yp1 - yq2 + yq1
      c              = xr1 - xq1
      d              = xr2 - xr1 - xq2 + xq1
      e              = yr1 - yq1
      f              = yr2 - yr1 - yq2 + yq1
      g              = xq1 - xp1
      h              = xq2 - xq1 - xp2 + xp1
      a'             = f*h + b*d
      b'             = f*g + e*h + a*d + b*c
      c'             = e*g + a*c
      d'             = b'*b' - 4 * a' * c'
  return (if d' < 0 then Nothing  else solveForCollision (gettime a b c d e f g h a' b' c' d') (time1, point1, segment1) (time2, point2, segment2) )

solveForCollision :: (Maybe Float, Maybe Float) -> (Float, Point 2 Float, LineSegment 2 () Float) -> (Float, Point 2 Float, LineSegment 2 () Float) -> Maybe (Float, Point 2 Float, Vector 2 Float)
solveForCollision (Nothing, Nothing) _ _ = Nothing

solveForCollision (Just t, Nothing) (time1, point1, segment1) (time2, point2, segment2) = let
      Point2 xp1 yp1 = point1
      Point2 xp2 yp2 = point2
      Point2 xq1 yq1 = segment1 ^. start . core
      Point2 xq2 yq2 = segment2 ^. start . core
      Point2 xr1 yr1 = segment1 ^. end . core
      Point2 xr2 yr2 = segment2 ^. end . core
      ptx             = xp1 + t * (xp2 - xp1)
      qtx             = xq1 + t * (xq2 - xq1)
      rtx             = xr1 + t * (xr2 - xr1)
      pty             = yp1 + t * (yp2 - yp1)
      qty             = yq1 + t * (yq2 - yq1)
      rty             = yr1 + t * (yr2 - yr1)
      f = if qtx == rtx then (pty - qty)/(rty - qty) else (ptx - qtx)/(rtx - qtx)
      segAngle = Prelude.atan ((qty - rty)/(qtx - rtx))
      reflectionMatrix = reflection segAngle
      xSegmentVelocity = f*(xr2 - xr1) + (1-f)*(xq2 - xq1)
      ySegmentVelocity = f*(yr2 - yr1) + (1-f)*(yq2 - yq1)
      xPointVelocity = xp2 - xp1
      yPointVelocity = yp2 - yp1
      xVelocityDifference = xPointVelocity - xSegmentVelocity
      yVelocityDifference = yPointVelocity - ySegmentVelocity
      reflectedVelocityDifference = transformBy reflectionMatrix  (Vector2 xVelocityDifference yVelocityDifference)
      Vector2 reflectedVelocityDifferencex reflectedVelocityDifferencey = reflectedVelocityDifference
      finalVelocity = Vector2 ((reflectedVelocityDifferencex + xSegmentVelocity)/(time2-time1)) ((reflectedVelocityDifferencey + ySegmentVelocity)/(time2-time1))
  in if f >= 0 && f <= 1 then Just (time1 + t * (time2 - time1), Point2 ptx pty, finalVelocity) else Nothing

solveForCollision (Just t1, Just t2) (time1, point1, segment1) (time2, point2, segment2)   = let
      Point2 xp1 yp1 = point1
      Point2 xp2 yp2 = point2
      Point2 xq1 yq1 = segment1 ^. start . core
      Point2 xq2 yq2 = segment2 ^. start . core
      Point2 xr1 yr1 = segment1 ^. end . core
      Point2 xr2 yr2 = segment2 ^. end . core
      ptx1             = xp1 + t1 * (xp2 - xp1)
      qtx1             = xq1 + t1 * (xq2 - xq1)
      rtx1             = xr1 + t1 * (xr2 - xr1)
      pty1             = yp1 + t1 * (yp2 - yp1)
      qty1             = yq1 + t1 * (yq2 - yq1)
      rty1             = yr1 + t1 * (yr2 - yr1)
      f1 = if qtx1 == rtx1 then (pty1 - qty1)/(rty1 - qty1) else (ptx1 - qtx1)/(rtx1 - qtx1)
      ptx2             = xp1 + t2 * (xp2 - xp1)
      qtx2             = xq1 + t2 * (xq2 - xq1)
      rtx2             = xr1 + t2 * (xr2 - xr1)
      pty2             = yp1 + t2 * (yp2 - yp1)
      qty2             = yq1 + t2 * (yq2 - yq1)
      rty2             = yr1 + t2 * (yr2 - yr1)
      f2 = if qtx2 == rtx2 then (pty2 - qty2)/(rty2 - qty2) else (ptx2 - qtx2)/(rtx2 - qtx2)
      f | (f1 >= 0 && f1 <= 1) = f1
        | (f2 >= 0 && f2 <= 1) = f2
        | otherwise = -1
      t | (f1 >= 0 && f1 <= 1) = t1
        | (f2 >= 0 && f2 <= 1) = t2
        | otherwise = -1
      segAngle
        | (f1 >= 0 && f1 <= 1) = Prelude.atan ((qty1 - rty1)/(qtx1 - rtx1))
        | (f2 >= 0 && f2 <= 1) = Prelude.atan ((qty2 - rty2)/(qtx2 - rtx2))
        | otherwise = -1
      reflectionMatrix = reflection segAngle
      xSegmentVelocity = f*(xr2 - xr1) + (1-f)*(xq2 - xq1)
      ySegmentVelocity = f*(yr2 - yr1) + (1-f)*(yq2 - yq1)
      xPointVelocity = xp2 - xp1
      yPointVelocity = yp2 - yp1
      xVelocityDifference = xPointVelocity - xSegmentVelocity
      yVelocityDifference = yPointVelocity - ySegmentVelocity
      reflectedVelocityDifference = transformBy reflectionMatrix  (Vector2 xVelocityDifference yVelocityDifference)
      Vector2 reflectedVelocityDifferencex reflectedVelocityDifferencey = reflectedVelocityDifference
      finalVelocity = Vector2 ((reflectedVelocityDifferencex + xSegmentVelocity)/(time2-time1)) ((reflectedVelocityDifferencey + ySegmentVelocity)/(time2-time1))
  in if f1 >= 0 && f1 <= 1 then Just (time1 + t1 * (time2 - time1), Point2 ptx1 pty1, finalVelocity) else (if f2 >= 0 && f2 <= 1 then Just (time1 + t2 * (time2 - time1), Point2 ptx2 pty2, Vector2 (-1) 0) else Nothing)


gettime :: Float  -> Float  -> Float  -> Float  -> Float  -> Float  -> Float  -> Float  -> Float  -> Float  -> Float -> Float -> (Maybe Float, Maybe Float)
gettime a b c d e f g h a' b' c' d' =
  if a' == 0 then (if -(a * c + e * g)/(f*g + e*h + a*d + c * b) >= 0 && -(a * c + e * g)/(f*g + e*h + a*d + c * b) <= 1 then (Just (-(a * c + e * g)/(f*g + e*h + a*d + c * b)), Nothing) else (Nothing, Nothing)) else
  let sol1 = (-b' - sqrt d') / (2 * a')
      sol2 = (-b' + sqrt d') / (2 * a')
  in (if sol1 <= 1 && sol1 >= 0 && sol2 <= 1 && sol2 >= 0 then (if sol1  <= sol2 then (Just sol1, Just sol2) else (Just sol2, Just sol1)) else (if sol1 <= 1 && sol1 >= 0  then (Just sol1, Nothing) else (if sol2 <= 1 && sol2 >= 0  then (Just sol2, Nothing) else (Nothing, Nothing))))


plan :: (Float, Arm) -- location and description of the arm
     -> (Point 2 Float, Vector 2 Float, Vector 2 Float) -- desired point of collision, orientation of bat, and velocity of the bat
     -> IO (Maybe ([Float], [Float])) -- output position and angular velocity of the arm

plan args1@(x, arm) args2@(pos, ori, vel) = ccdprecompute args1 args2

type Angle = Float -- angle (should be in radians)
type AngulVel = Float --angular speed (should be in radians per second)
type CustomMatrix = [[Float]] -- custom matrix type for calculations, inner list is a row

ccdprecompute :: (Float, Arm) -- location and description of the arm
     -> (Point 2 Float, Vector 2 Float, Vector 2 Float) -- desired point of collision, orientation of bat, and velocity of the bat
     -> IO (Maybe ([Angle], [AngulVel])) -- output position and angular velocity of the arm

--precompute the two points at which the bat could end and use that as target for ccd for rest of the arm

ccdprecompute (x, arm) (pos, ori, vel@(Vector2 velx vely)) = do
  let
    Link _ batlength : _ : restofarm = reverse arm
    Point2 posx posy = pos
    Vector2 orix oriy = signorm ori
    batEnd1 = Point2 (posx - (oriy * batlength / 2)) (posy + (orix * batlength / 2))
    batEnd2 = Point2 (posx + (oriy * batlength / 2)) (posy - (orix * batlength / 2))
    batAngle1 = atan2 (-oriy) (-orix)
    batAngle2 = atan2 oriy orix
    armangles1 = ccdentry (x, restofarm) batEnd1
    armangles2 = ccdentry (x, restofarm) batEnd2
    armAnglesResult = reverse $ returnFinalAnswer  (armangles1, batAngle1) (armangles2, batAngle2)
    vectorizedresultarm = halfFirstVector $ vectorizeArm $ reverse $ evaluateU (combineAnglesInArm arm armAnglesResult)
    perpendiculars = map (\(Vector2 x y)-> Vector2 (-y) x) $ reverse $ computePerpendiculars (Prelude.init vectorizedresultarm)
    perpendicularsx = map (\(Vector2 x _) -> x) perpendiculars
    perpendicularsy = map (\(Vector2 _ y) -> y) perpendiculars
    jacobian = [perpendicularsx, perpendicularsy]
    jacobianInverse = if length (Prelude.head jacobian) == 2
      then inverseMatrix22 jacobian
      else pseudoInverseMatrix jacobian
    wantedAngulVelUnformatted = multiplyMatrix jacobianInverse [[velx],[vely]]
    wantedAngulVelFormatted = Prelude.head $ transposeMatrix wantedAngulVelUnformatted
  return (if not $ null armAnglesResult then Just (armAnglesResult, wantedAngulVelFormatted) else Nothing)

pseudoInverseMatrix :: CustomMatrix -> CustomMatrix
pseudoInverseMatrix matrix = multiplyMatrix
  (transposeMatrix matrix)
  (inverseMatrix22 (multiplyMatrix matrix (transposeMatrix matrix)))

transposeMatrix:: CustomMatrix -> CustomMatrix --https://stackoverflow.com/questions/2578930/understanding-this-matrix-transposition-function-in-haskell/2578979
transposeMatrix ([]:_) = []
transposeMatrix matrix = map Prelude.head matrix : transposeMatrix (map tail matrix)

multiplyMatrix :: CustomMatrix -> CustomMatrix -> CustomMatrix
multiplyMatrix a b = [[ sum $ zipWith (*) ar bc | bc <- transposeMatrix b] | ar <- a] --https://samprocter.com/2012/11/matrix-multiplication-in-haskell/

inverseMatrix22 :: CustomMatrix -> CustomMatrix
inverseMatrix22 [row1@[a, b], row2@[c, d]] = let determinant = a*d - b*c
  in map (map (/determinant)) [[d, -b], [-c, a]]
inverseMatrix22 matrix = trace ("This wasn't a 2x2 matrix" ++ show matrix) []

halfFirstVector :: [Vector 2 Float] -> [Vector 2 Float]
halfFirstVector (x : xs) = (x ^/ 2) : xs
halfFirstVector _ = []

computePerpendiculars :: [Vector 2 Float] -> [Vector 2 Float]
computePerpendiculars (curvec : nextvec : rest) = curvec : computePerpendiculars (curvec ^+^ nextvec : rest)
computePerpendiculars [curvec] = [curvec]
computePerpendiculars [] = []

combineAnglesInArm :: Arm -> [Angle] -> Arm
combineAnglesInArm (Joint c _ : restarm) (a : rest) = Joint c a : combineAnglesInArm restarm rest
combineAnglesInArm (link@(Link _ _) : restarm) angles = link : combineAnglesInArm restarm angles
combineAnglesInArm _ _ = []

ccdentry :: (Float, Arm) -- location and description of the arm
     -> Point 2 Float -- ccd target point
     -> Maybe [Angle] -- output position of the arm

ccdentry (x, arm) target@(Point2 tx ty) = let resultarm = ccd 100 (x, arm) target
  in let result = armToAngles resultarm
  in let (Point2 px py) : _ = reverse (evaluateU (reverse resultarm))
  in let distance = sqrt((square (tx - px - x)) + (square (ty - py)))
  in if distance < 0.0001 then Just result else Nothing --make check correct

returnFinalAnswer :: (Maybe [Float], Angle) -> (Maybe [Float], Angle) -> [Float]
returnFinalAnswer (Just x, a) _ =  modAngle (a - sum x) : x
returnFinalAnswer  _ (Just x, a) = modAngle (a - sum x) : x
returnFinalAnswer _ _ = []

square :: Float -> Float
square x = x*x

armToAngles :: Arm -> [Angle]
armToAngles [] = []
armToAngles (Joint _ a : x) = a : armToAngles x
armToAngles (Link _ _ : x) = armToAngles x

ccd :: Int
     -> (Float, Arm) -- location and description of the arm
     -> Point 2 Float -- ccd target point
     -> Arm -- output position of the arm

ccd iter (x, arm) target@(Point2 px py) =
  let pos = Point2 (px - x) py
  in let evalArm = reverse (evaluateU (reverse arm))
  in let vectArm = vectorizeArm evalArm
  in let iterresult = ccditer vectArm (drop 1 evalArm) pos arm
  in (if iter <= 1 then iterresult else ccd (iter-1) (x, iterresult) target)


fixEvaluatedArm :: [a] -> [a]
fixEvaluatedArm (x : _ : rest) = x : fixEvaluatedArm rest
fixEvaluatedArm _ = []

applyMaybeOrReturnNot :: (a -> Maybe b) -> Maybe a -> Maybe b
applyMaybeOrReturnNot func (Just x) = func x
applyMaybeOrReturnNot _ _ = Nothing

vectorizeArm :: [Point 2 Float] -> [Vector 2 Float]
vectorizeArm ((Point2 x y) : (Point2 xs ys) : rest) = Vector2 (x - xs) (y - ys) : vectorizeArm (Point2 xs ys: rest)
vectorizeArm [_] = []
vectorizeArm [] = []

--iteration function that will call itself untill iter float is 0
ccditer :: [Vector 2 Float] -- location and length of the arm links
     ->[Point 2 Float]
     -> Point 2 Float -- ccd target point
     -> Arm -- intermediate position of the arm
     -> Arm -- output position of the arm

ccditer _ _ _ [] = []
ccditer _ _ _ [Link a b] = [Link a b]
ccditer (armend@(Vector2 armx army) : nextarm : armrest) (hingepos@(Point2 hingex hingey) : hingerest) pos@(Point2 x y) (link@(Link _ _) : (Joint c curangle) : restcurarm) = let
  wantedDirVec@(Vector2 vecx vecy) = Vector2 (x - hingex) (y - hingey)
  in let wantedRotation = atan2  vecy vecx - atan2 army armx
  in let newAngle = modAngle (curangle + wantedRotation)
  in let newArmend = (sqrt (armx*armx + army*army) *^ signorm wantedDirVec) ^+^ nextarm
  in link : Joint c newAngle : ccditer (newArmend : armrest) hingerest pos restcurarm

ccditer _ _ _ _ = []