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
  for (int i; i < num_particles; i++){
    std::cout << "weight[i] is " << i << std::endl;
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
    std::cout << "weight[i] is " << weights[i]<< std::endl;
  }
  is_initialized = true;
  std::cout << "init is done " << std::endl;
  
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

  for(unsigned int i = 0; i < particles.size(); i++){
    
    //when yaw_rate is close to zero, use linear model
    if(fabs(yaw_rate) < 0.00001){
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      //particles[i].theta does not change
      std::cout << "yaw_rate is close to 0"<< std::endl;
    }
    
    //when yaw_rate != zero, use general model
    else{
      particles[i].x += velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;
      std::cout << "yaw_rate is normal"<<std::endl;
    }
    //prediction, add gaussian noise for x
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    particles[i].x = dist_x(gen);
    //prediction, add gaussian noise for y
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    particles[i].y = dist_y(gen);
    //prediction, add gaussian noise for theta
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    particles[i].theta = dist_theta(gen);
  }
std::cout << "prediction is done " << std::endl;
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
  //predicted is a set of landmarks within sensor range
  for(unsigned int i = 0; i <predicted.size(); i++){
    //initialize the distance with infinity 
    double distance = std::numeric_limits<double>::infinity();
    int closest = 0;
    for(unsigned int j = 0; j<observations.size(); j++){
      //loop through to find shortest distance with the helper function dist()
      double dist_ij = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
      if(dist_ij < distance){
        distance = dist_ij;
        closest = j;
      }
    }
    predicted[i].id = observations[closest].id;
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
  //x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
  //y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
  double weight_sum = 0.0;
  for(unsigned int i = 0; i < particles.size(); i++){
    double weight = 1.0;
    
    //vector<LandmarkObs> observations_map;
    for(unsigned int j = 0; j< observations.size(); j++){
      LandmarkObs observation_map;
      observation_map.x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      observation_map.y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
	  observation_map.id = observations[j].id;
      //observations_map.push_back(observation_map);
      
      int closest = 0;
      for(unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++){
        double particle_landmark_distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f); //check if within sensor range
        double distance = std::numeric_limits<double>::infinity();//init the distance between observation and landmark with infinity number
        std::cout<<"distance between observation and landmark is " << distance << std::endl; //debug purpose

        // Step2. ENSURE map landmarks are inside sensor range
        if(particle_landmark_distance <= sensor_range){
          // Step3. Nearest Neighbor Data Association
          double landmark_obs_distance = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, observation_map.x, observation_map.y);
          if(landmark_obs_distance < distance ){
            closest = k;
            distance = landmark_obs_distance;
          }
        }
      }
      observation_map.id = map_landmarks.landmark_list[closest].id_i;
      weight *=  multiv_prob(std_landmark[0], std_landmark[1], observation_map.x, observation_map.y, map_landmarks.landmark_list[closest].x_f, map_landmarks.landmark_list[closest].y_f);
    }
    //update the weight
    particles[i].weight = weight;
    weights[i] = particles[i].weight;
    weight_sum += weight;
    std::cout << "weight[i] is " << weights[i] << std::endl;
    // Step2. ENSURE map landmarks are inside sensor range
    /*vector<LandmarkObs> predicted;
    for(unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++){
      double particle_landmark_distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f); 
      LandmarkObs landmark;
      
      if(particle_landmark_distance <= sensor_range){
        landmark.id = map_landmarks.landmark_list[k].id_i;
        landmark.x = map_landmarks.landmark_list[k].x_f;
        landmark.y = map_landmarks.landmark_list[k].y_f;
        predicted.push_back(landmark);
     } 
    }*/
      // Step3. Nearest Neighbor Data Association
      //Here dataAssociation function needs to be called.
      //dataAssociation(predicted, observations_map);
    
      // Step4. Compute WEIGHT of particle
      //double weight = 1;
     /* for(unsigned int q; q <predicted.size(); q++){
        for(unsigned int w; w< observations_map.size(); w++){
          if(predicted[q].id == observations_map[w].id){
            weight *=  multiv_prob(std_landmark[0], std_landmark[1], observations_map[w].x, observations_map[w].y, predicted[q].x, predicted[q].y);
            
          }
        }
      }
      particles[i].weight = weight;
      weights[i] = particles[i].weight;*/
  }
  if(fabs(weight_sum >0.0)){
    for(unsigned int i = 0; i < weights.size(); i++){
      weights[i] = weights[i]/weight_sum;
    }
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