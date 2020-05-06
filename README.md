# Autonomous-IterativeClosePoint

Helper files have been added in the helper files. ROS package modified 

What is -IterativeClosePoint


One of the most fundamental algorithms of
localization, scan matching, is implemented. It uses the Iterative Closest Point algorithm.
Concepts Involved in this assignment

1. Localization
2. Odometry Estimation
3. Convex optimization
4. C++ OOP
5. Quadratic Programming


ICP Algorithm


Potentially sample Points
1. Determine corresponding points
2 Potentially weight / reject pairs
3. Compute rotation R, translation t (e.g. SVD)
4. Apply R and t to all points of the set to be
registered
5. Compute the error E(R,t)
6. If error decreased and error > threshold
7.  Repeat to determine correspondences etc.
8. Stop and output final alignment, otherwise 



Data Structures used within the files 

The skeleton uses a few data structures which are dened in the dierent program les. You are free to use the
data structures provided or make new ones of your own. The functions however should be implemented as dened.
1. Point
The struct the radial distance and angular distance of a point from the car frame. In this struct there are few
functions implemented which can be used to derive other information from the points:
    distToPoint(point P) : find the distance to another point
    distToPoint2(point P) : find the square of the distance to another point
    radialGap(point P) : find the radial gap to another point
    getx() : get the x coordinate of the point
    gety() : get the y coordinate of the point
    wrapTheta() : wrap theta around 2pi during rotation
    rotate(phi) : rotate the point by angle phi
    translate(x,y) : translate a point by distance x and y
    getVector(): get the vector (in x and y)


2. Correspondence
This struct stores the correspondence values which you nd through the fast search algorithm. It contains
the transformed points: P, original points : Po , rst best point: Pj1, second best point pj2.
    getNormal() Get normal of the correspondence
    getPiGeo() get correspondence point as a geometry message
    getPiVec() get correspondence point as a vector message


3. Transform
The struct stores the transform which you calculate through optimization at every step. it contains the x
translation, y translation and the theta rotation.
    apply(point P): apply the transform to a provided point
    getMatrix() : get the transformation matrix from the found transform


