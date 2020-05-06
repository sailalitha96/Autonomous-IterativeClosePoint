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

