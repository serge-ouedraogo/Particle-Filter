/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using namespace std;




void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   
   num_particles = 100;  // TODO: Set the number of particles
   random_device rd;
   default_random_engine gen(rd());
   
   weights.resize(num_particles);
   particles.resize(num_particles);
 
   normal_distribution<double> x_distribution(x, std[0]);
   normal_distribution<double> y_distribution(y, std[1]);
   normal_distribution<double> theta_distribution(theta, std[2]);
   
  // Initialized particles coordinates with randomly generated numbers;
  for(int i = 0; i < num_particles; ++i)
  {
    particles[i].id = i;
    particles[i].x = x_distribution(gen);
    particles[i].y = y_distribution(gen);
    particles[i].theta = theta_distribution(gen);
    particles[i].weight = 1.0;
  }
   
  is_initialized = true;
  
  }

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // Implement motion model on each particle
  
  
  default_random_engine gen;
 
  normal_distribution<double> x_distribution(0, std_pos[0]);
  normal_distribution<double> y_distribution(0, std_pos[1]);
  normal_distribution<double> theta_distribution(0, std_pos[2]);
  
  for(int i = 0; i < num_particles; ++ i)
  {
    if(fabs(yaw_rate) < 0.00001)
    {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t *  sin(particles[i].theta);
    }
    else
    {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));     
      particles[i].y +=  velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
     //add error terms
     particles[i].x += x_distribution(gen);
     particles[i].y += y_distribution(gen);
     particles[i].theta += theta_distribution(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) 
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
    //loop through all the observed landmarks
   
    //cout << "Observation SIZE: " << observations.size() << endl; 
    
    for(unsigned i =0; i < observations.size(); ++i)
    {   
      int nearest_particle_id = -10; // arbitrary initialised the nearest particle. 
      double min_dist = 1E9;
      // For each observed landmark loop through all the possible (predicted) landmarks
      for(unsigned j =0; j < predicted.size(); ++j)
      {        
        // identify the landmarks between the two groups (predicted and observed ) closest to each other
        double lm_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
        if(lm_dist < min_dist)
        {
          min_dist = lm_dist;
          nearest_particle_id = j; 
        }
      }
     
      observations[i].id = nearest_particle_id;
      //cout << "Nearest Particle: " << observations[i].id << endl; 
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
   // A sensor was used to take a measurement so let's update our knowledge.
  
  auto norm_factor = 2*M_PI* std_landmark[0] *  std_landmark[1]; 
  auto cov_x =  std_landmark[0]  * std_landmark[0];
  auto cov_y =  std_landmark[1] *  std_landmark[1];
  
 
  for(int i = 0; i < num_particles; ++i)
  {
    vector<LandmarkObs> predicted;
    for(unsigned j = 0; j < map_landmarks.landmark_list.size(); ++j)
    {
      double lmdist = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y);
      
      // Ensure the sensor range covers all the predicted landmarks.
      if(lmdist <= sensor_range)
      {
        predicted.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
      }
    }
    
    if(predicted.size()==0)
    {
      particles[i].weight = 0; 
      weights[i] = 0;
    }
    else
    {
      //cout << "Predicted SIZE: " << predicted.size() << endl; 
      // Transform coordinates of observed landmarks from Car ref. frame to map ref frame
      vector<LandmarkObs> trans_obs;
      for(unsigned p = 0; p < observations.size(); ++p)
      {
        double tx = cos(particles[i].theta) * observations[p].x - sin(particles[i].theta) * observations[p].y + particles[i].x;
        double ty = sin(particles[i].theta) * observations[p].x + cos(particles[i].theta) * observations[p].y + particles[i].y;
        trans_obs.push_back(LandmarkObs{observations[p].id, tx, ty});
      }
      
      //cout << "Trans_OBS SIZE: " << trans_obs.size() << endl; 
      dataAssociation(predicted, trans_obs);
      double prob = 1.0;
      
      for( unsigned q = 0; q < trans_obs.size(); ++q)
      {
        auto dx = trans_obs[q].x - predicted[trans_obs[q].id].x;
        auto dy = trans_obs[q].y - predicted[trans_obs[q].id].y;
        
        prob *= exp(-(dx * dx / (2 * cov_x) + dy * dy / (2 * cov_y))) / norm_factor;
        
       
      }
      //cout << "prob: " << prob << endl;  
      particles[i].weight = prob;
      weights[i] = prob;
    }
  }
}

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  random_device rd;
  default_random_engine gen(rd());
  discrete_distribution<int> p_filter(weights.begin(), weights.end());
  vector<Particle> newparticles(num_particles);
  for(int i = 0; i < num_particles; ++i)
  {
    newparticles[i] = particles[p_filter(gen)];
  }
  
  //Use move to replace old particles with new particles
  particles = move(newparticles); 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}