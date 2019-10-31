#include "scan_matching_skeleton/correspond.h"
#include "cmath"
#include "ros/ros.h"

using namespace std;

const int UP_SMALL = 0;
const int UP_BIG = 1;
const int DOWN_SMALL = 2;
const int DOWN_BIG = 3;

void getNaiveCorrespondence(vector<Point>& old_points, vector<Point>& trans_points, vector<Point>& points,
                        vector< vector<int> >& jump_table, vector<Correspondence>& c, float prob){

      c.clear();
      int last_best = -1;
      const int n = trans_points.size();
      const int m = old_points.size();
      float min_dist = 100000.00;
      int min_index = 0;
      int second_min_index = 0;

      //Do for each point
      for(int i = 0; i<n; ++i){
      	min_dist = 100000.00;
        for(int j = 0; j<m; ++j){
          float dist = old_points[j].distToPoint2(&trans_points[i]);
          if(dist<min_dist){
            min_dist = dist;
            min_index = j;

  			if(min_index==0) { second_min_index = min_index+1;} 
  			else {second_min_index = min_index-1;}
          }
        }
        c.push_back(Correspondence(&trans_points[i], &points[i], &old_points[min_index], &old_points[second_min_index]));
      }
}

void getCorrespondence(vector<Point>& old_points, vector<Point>& trans_points, vector<Point>& points,
                        vector< vector<int> >& jump_table, vector<Correspondence>& c, float prob){

  // Written with inspiration from: https://github.com/AndreaCensi/gpc/blob/master/c/gpc.c
  // use helper functions and structs in transform.h and correspond.h
  // input : old_points : vector of struct points containing the old points (points of the previous frame)
  // input : trans_points : vector of struct points containing the new points transformed to the previous frame using the current estimated transform
  // input : points : vector of struct points containing the new points
  // input : jump_table : jump table computed using the helper functions from the transformed and old points
  // input : c: vector of struct correspondences . This is a refernece which needs to be updated in place and return the new correspondences to calculate the transforms.
  // output : c; update the correspondence vector in place which is provided as a reference. you need to find the index of the best and the second best point. 
  //Initializecorrespondences

  c.clear();
  int last_best = -1;
  const int n = trans_points.size();
  const int m = old_points.size();
  // iterate each point 
  for(int i = 0; i<n; ++i){
    int begin_at=0;
    // int start_at = (last_best<0)?(last_best+1):last_best;
    if(last_best<0)
    {
      begin_at = last_best+1;
    }
    else{
      begin_at = last_best;
    }

	   int min_idx = begin_at;
	   int sec_min_idx;
	   float min_dist = 10000.0;

    bool upward= false;
    bool downward = false;

    // start from last best, search upwards
    int curr = begin_at;

    while(!upward && curr<m){
    	float dt = trans_points[i].distToPoint2(&old_points[curr]);  

    	if (dt<min_dist){
    		  min_dist = dt;// set dist to curr cal 
    		  min_idx = curr;// set idx , you didnt do that (fix it )

    		if(min_idx==0) { 
            sec_min_idx = min_idx+1;
          } 
  			else { 
          sec_min_idx= min_idx-1;
        }
    	}
      // declare the dfference in the angle for up direcion 
    	float angle_dif= old_points[curr].theta - trans_points[i].theta;

    	if (trans_points[i].r * sin(angle_dif) > min_dist)
      { 
    		    upward = true;// set the contratry 
    	}
    	if(old_points[curr].r < trans_points[i].r)
      {
    		curr = jump_table[curr][UP_BIG];// jump advanced in the table 
    	}
    	else{
    		curr = jump_table[curr][UP_SMALL];// jump advanced on the tabl e
    	}
    }
    curr = begin_at;

    while(!downward && curr>=0){
    	float dt = trans_points[i].distToPoint2(&old_points[curr]);    	
    	if (dt<min_dist){
    		min_dist = dt;
    		min_idx = curr;

    		if(min_idx==0) { sec_min_idx= min_idx+1;} 
  			else { sec_min_idx= min_idx-1;}
    	}
    	float angle_diff1 = trans_points[i].theta - old_points[curr].theta;

    	if (trans_points[i].r * sin(angle_diff1) > min_dist){ 
    		downward = true;
    	}
    	// advance based on jump table
    	if(old_points[curr].r < trans_points[i].r){
    		curr = jump_table[curr][DOWN_BIG];
    	}
    	else{
    		curr = jump_table[curr][DOWN_SMALL];
    	}
    }

    c.push_back(Correspondence(&trans_points[i], &points[i], &old_points[min_idx], &old_points[sec_min_idx]));
    last_best = min_idx;
    }
  }


void computeJump(vector< vector<int> >& table, vector<Point>& points){
  table.clear();
  int n = points.size();
  for(int i = 0; i<n; ++i){
    vector<int> v = {n,n,-1,-1};
    for(int j = i+1; j<n; ++j){
      if(points[j].r<points[i].r){
        v[UP_SMALL] = j;
        break;
      }
    }

    for(int j = i+1; j<n-1; ++j){
      if(points[j].r>points[i].r){
        v[UP_BIG] = j;
        break;
      }
    }
    for(int j = i-1; j>=0; --j){
      if(points[j].r<points[i].r){
        v[DOWN_SMALL] = j;
        break;
      }
    }
    for(int j = i-1; j>=0; --j){
      if(points[j].r>points[i].r){
        v[DOWN_BIG] = j;
        break;
      }
    }
    table.push_back(v);
  }
}
