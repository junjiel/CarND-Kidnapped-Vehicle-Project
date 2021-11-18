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
#include <limits>

#include "helper_functions.h"
using std::string;
using std::vector;
using std::normal_distribution;

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
  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  for (int i = 0; i < num_particles; i++){
    Particle p;
    p.id = i;
    // Sample from these normal distributions like this: 
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1;
    particles.push_back(p);
    weights.push_back(p.weight);
  }
  is_initialized = true;
  //std::cout << "init is done " << std::endl;
  
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
  std::default_random_engine gen;
  //create gaussian distribution with mean 0
  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);
  for(unsigned int i = 0; i < particles.size(); i++){
    
    //when yaw_rate is close to zero, use linear model
    if(fabs(yaw_rate) < 0.00001){
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      //particles[i].theta does not change
      //std::cout << "yaw_rate is close to 0"<< std::endl;
    }
    
    //when yaw_rate != zero, use general model
    else{
      particles[i].x += velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;
      //std::cout << "yaw_rate is normal"<<std::endl;
    }
    //prediction, add gaussian noise for x
    particles[i].x += dist_x(gen);
    //prediction, add gaussian noise for y
    particles[i].y += dist_y(gen);
    //prediction, add gaussian noise for theta
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for(unsigned int i = 0; i <observations.size(); i++){
    //initialize the distance with max 
    double distance = std::numeric_limits<double>::max();
    int closest = -1;
    for(unsigned int j = 0; j<predicted.size(); j++){
      //loop through to find shortest distance with the helper function dist()
      double dist_ij = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      if(dist_ij < distance){
        distance = dist_ij;
        closest = predicted[j].id;
      }
    }
    observations[i].id = closest;
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
  
  // Step1. TRANSFORM each observation marker from the vehicle's coordinates to the map's coordinates
  //for each particle, it should transform the observations to map coordinates, so this transformation should be within for loop of particles
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
  
  for(unsigned int i = 0; i < particles.size(); i++){
    // Step1. TRANSFORM each observation marker from the vehicle's coordinates to the map's coordinates
    vector<LandmarkObs> observations_map;
    for(unsigned int j = 0; j< observations.size(); j++){
      double x_map = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      double y_map = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
      observations_map.push_back (LandmarkObs{observations[j].id, x_map, y_map});
    }
    
    // Step2. ENSURE map landmarks are inside sensor range
    vector<LandmarkObs> predicted;
    for(unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++){
      int land_ID = map_landmarks.landmark_list[k].id_i;
      float land_x = map_landmarks.landmark_list[k].x_f;
      float land_y = map_landmarks.landmark_list[k].y_f;
      double particle_landmark_distance = dist(particles[i].x, particles[i].y, land_x, land_y); 
      
      if(particle_landmark_distance <= sensor_range){
        predicted.push_back(LandmarkObs{land_ID, land_x, land_y});
     } 
    }
    
    // Step3. Nearest Neighbor Data Association
    //Here dataAssociation function needs to be called.
    dataAssociation(predicted, observations_map);
    
    // Step4. Compute WEIGHT of particle
    vector<int> association;
    vector<double> sense_x;
	vector<double> sense_y;
    particles[i].weight = 1.0;
    //double p_weight = 1.0;
    double map_x, map_y, mu_x, mu_y;
   for (unsigned int t = 0; t < observations_map.size(); ++t)
    {
      map_x =  observations_map[t].x;
      map_y =  observations_map[t].y;
      for (unsigned int p = 0; p < predicted.size(); ++p)
      {
        // Associate prediction with transformed observation
        if (predicted[p].id == observations_map[t].id)
        {
          mu_x = predicted[p].x;
          mu_y = predicted[p].y;
          //p_weight = multiv_prob(std_landmark[0], std_landmark[1], map_x, map_y, mu_x, mu_y);
        }
      }
      // Compute exponent
      double exponent = (0.5*pow( (map_x - mu_x)/sig_x, 2.0 )+0.5*pow( (map_y - mu_y)/sig_y, 2.0 ));
      // Compute weight using normalization terms and exponent
      double p_weight = gauss_norm * exp(- exponent);
      
      particles[i].weight *= p_weight;
      
      // Append particle associations
      association.push_back(observations_map[t].id);
      sense_x.push_back(observations_map[t].x);
      sense_y.push_back(observations_map[t].y);
      }
    weights[i] = particles[i].weight;
    // For blue lasers (belongs to the best particle )
    SetAssociations(particles[i], association, sense_x, sense_y);    
 /*
    for(unsigned int w; w< observations_map.size(); w++){
      for(unsigned int q; q <predicted.size(); q++){
        if(predicted[q].id == observations_map[w].id){
          weight =  multiv_prob(std_landmark[0], std_landmark[1], observations_map[w].x, observations_map[w].y, predicted[q].x, predicted[q].y);          
          }
        }
      particles[i].weight *= weight;*/
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
  std::default_random_engine gen;
  // Create uniform distribution
  std::discrete_distribution<> dist_weighted(weights.begin(), weights.end());
  vector<Particle> particles_sampled;
    
  for (unsigned int i = 0; i < particles.size() ; ++i) {
    int sample_index = dist_weighted(gen);
    particles_sampled.push_back(particles[sample_index]);
  }
  particles = particles_sampled;  

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
  s = s.substr(0, s.length()-1);  // get ridö of the trailing space
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