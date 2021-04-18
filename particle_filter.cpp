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
#include <chrono>
#include <cmath>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    // This line creates a normal (Gaussian) distribution for x
  cout<<"Initializing with the following input: "<<"X: "<<x<<" Y: "<<y<<" Theta: "<<theta<<endl;
  
  std::default_random_engine gen;
  std::normal_distribution<double> distX(x, std[0]);
  std::normal_distribution<double> distY(y, std[1]);
  std::normal_distribution<double> distT(theta, std[2]);
  num_particles = 10;  // TODO: Set the number of particles
  for(int i = 0; i<num_particles; ++i){
    Particle P;
    P.id = i;
    P.x = distX(gen);
    P.y = distY(gen);
    P.theta = distT(gen);
    P.weight = 1.0;
    P.associations = {};
    P.sense_x = {};
    P.sense_y = {};
    particles.push_back(P);
    weights.push_back(P.weight);
    cout<<"Initializing particle "<<i<<": ("<<P.x<<", "<<P.y<<") @ "<<P.theta<<endl;
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
  std::default_random_engine gen;

  for(int i = 0; i<num_particles; ++i){
    //cout<<"Predicting particle "<<i<<": ("<<particles[i].x<<", "<<particles[i].y<<") @"<<particles[i].theta<<" --> ";
    std::normal_distribution<double> distX(0, std_pos[0]);
    std::normal_distribution<double> distY(0, std_pos[1]);
    std::normal_distribution<double> distT(0, std_pos[2]);
    if(yaw_rate != 0){
      
      particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + 0.01*distX(gen);
      particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + 0.01*distY(gen);
      particles[i].theta += yaw_rate*delta_t + 0.01*distT(gen);
    }
    else{
      particles[i].x += velocity*delta_t*cos(particles[i].theta) + 0.01*distX(gen);
      particles[i].y += velocity*delta_t*sin(particles[i].theta) + 0.01*distY(gen);
      particles[i].theta += 0.01*distT(gen);

    }
    //cout<<"("<<particles[i].x<<", "<<particles[i].y<<") @"<<particles[i].theta<<endl;
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

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {

	double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double totalWeight = 0;

    for(int i=0; i<num_particles; ++i){
        Particle *p = &particles[i];
        double W = 1.0;

        // convert observation from vehicle's to map's coordinate system
        for(int j=0; j<observations.size(); ++j){
            LandmarkObs currentObs = observations[j];
            LandmarkObs transformedObs;

            transformedObs.x = (currentObs.x * cos(p->theta)) - (currentObs.y * sin(p->theta)) + p->x;
            transformedObs.y = (currentObs.x * sin(p->theta)) + (currentObs.y * cos(p->theta)) + p->y;
            transformedObs.id = currentObs.id;

            // find the predicted measurement that is closest to each observed measurement and assign
            // the observed measurement to this particular landmark
            Map::single_landmark_s ldmk;
            double distance_min = std::numeric_limits<double>::max();

            for(int k=0; k<map_landmarks.landmark_list.size(); ++k){
                Map::single_landmark_s mapLandmark = map_landmarks.landmark_list[k];
                double distance = dist(transformedObs.x, transformedObs.y, mapLandmark.x_f, mapLandmark.y_f);
                if(distance < distance_min){
                    distance_min = distance;
                    ldmk = mapLandmark;
                }
            }

            // update weights using Multivariate Gaussian Distribution
            // equation given in Transformations and Associations Quiz
            double num = exp(-0.5 * (pow((transformedObs.x - ldmk.x_f), 2) / pow(std_x, 2) + pow((transformedObs.y - ldmk.y_f), 2) / pow(std_y, 2)));
            double denom = 2 * M_PI * std_x * std_y;
            W *= num/denom;
        }
        totalWeight += W;
        p->weight = W;
    }
    // normalize weights to bring them in (0, 1]
    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];
        p->weight /= totalWeight;
        weights[i] = p->weight;
    }
}





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
  /*
  double totalWeight = 0.0;
  cout<<"UPDATING WEIGHTS"<<endl;
  for(int i = 0; i<num_particles; ++i){
    //cout<<"Particle "<<i<<" ("<<particles[i].x<<", "<<particles[i].y<<") @"<<particles[i].theta<<endl;
    // Transform measured observations from vehicle frame to map frame
   
    //vector<LandmarkObs> transformedObs;
    double wt = 1.0;
    for(unsigned int j = 0; j<observations.size(); ++j){
      LandmarkObs newObs;
      newObs.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[i].y*sin(particles[i].theta);
      newObs.y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[i].y*cos(particles[i].theta);
      newObs.id = observations[j].id;
      //transformedObs.push_back(newObs);
    
      // Grab list of landmarks within sensor range
      vector<LandmarkObs> observable;
      for(unsigned int j = 0; j<map_landmarks.landmark_list.size(); ++j){
        LandmarkObs obs;
        double dx = map_landmarks.landmark_list[j].x_f - particles[i].x;
        double dy = map_landmarks.landmark_list[j].y_f - particles[i].y;
        //cout<<"\tDistance to Landmark "<<j<<": "<<sqrt(dx*dx + dy*dy)<<endl;
        if( sqrt(dx*dx + dy*dy) <= sensor_range ){
          obs.x = map_landmarks.landmark_list[j].x_f;
          obs.y = map_landmarks.landmark_list[j].y_f;
          obs.id = map_landmarks.landmark_list[j].id_i;
          observable.push_back(obs);
        }
      }

      //cout<<"# of Observable targets: "<<observable.size()<<endl;
      // Associate each transformed observation with an observable landmark
      Map::single_landmark_s ldmk;
      double association_distance = std::numeric_limits<double>::max();
      for(int j = 0; j<observable.size(); ++j){
        Map::single_landmark_s LDMK = map_landmarks.landmark_list[j];
        double testDistance = dist(newObs.x, newObs.y, LDMK.x_f, LDMK.y_f);
        if(testDistance < association_distance){
          association_distance = testDistance;
          ldmk = LDMK;
        }
      }
      
      double numer = exp(-0.5* (pow((newObs.x - ldmk.x_f), 2) / pow(std_landmark[0], 2) + pow((newObs.y - ldmk.y_f), 2) / pow(std_landmark[1], 2)));
      double denom = 2 * M_PI * std_landmark[0] * std_landmark[1];
      wt *= numer/denom;      
    }
    totalWeight += wt;
    particles[i].weight = wt;
  }
  
  cout<<"TOTAL WEIGHT: "<<totalWeight<<endl;
  
  // Normalize the weights of all particles
  //cout<<"NORMALIZING WEIGHTS"<<endl;
  cout<<"============================ Updated Weights ======================================"<<endl;
  for(int i = 0; i<num_particles; ++i){
    particles[i].weight /= totalWeight;
    cout<<"Particle "<<i<<" ("<<particles[i].x<<", "<<particles[i].y<<") @"<<particles[i].theta<<"\tW: "<<particles[i].weight<<"\tTotal Weight: "<<totalWeight<<endl;
    //if(std::isinf(particles[i].weight)){ particles[i].weight = 1; }
  }
}
*/
  
  



void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  std::normal_distribution<double> distX(0, 0.3);
  std::default_random_engine gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  vector<Particle> newParticles;
  //cout<<"Testing random generator: "<<gen()<<endl;
  //cout<<"RAND_MAX: "<<RAND_MAX<<endl;
  double rr = (double) gen()/RAND_MAX;
  //cout<<"rr: "<<rr<<" | "<<(double) (gen()/RAND_MAX)<<endl;
  int index = (int) rr*num_particles;
  //double maxWeight = max_element(weights);
  double maxWeight = 0.0;
  for(int i = 0; i<num_particles; ++i){if(particles[i].weight>maxWeight){maxWeight = particles[i].weight;}}
  
  double beta = 0.0;
  for(int i = 0; i<num_particles; ++i){
    rr = (double) gen()/RAND_MAX;
    //cout<<"2  rr: "<<rr<<endl;
    beta += rr*2.0*maxWeight;
    //cout<<"BETA: "<< beta;
    while(beta > particles[index].weight){
      //cout<<" INDEX: "<<index<<" | BETA: "<<beta<<endl;
      beta -= particles[index].weight;
      index = (index + 1) % num_particles;
      //cout<<" | "<<beta<<", "<<index;
    }
    //cout<<"\tResampling #"<<index<<endl;
    Particle newParticle = particles[index];
    newParticle.x += 0.1*distX(gen);
    newParticle.y += 0.1*distX(gen);
    newParticles.push_back(newParticle);
  }
  particles = newParticles;
  cout<<"=================== Resampled particles ==========================="<<endl;
  for(int i = 0; i<num_particles; ++i){
    cout<<"Particle "<<i<<": ("<<particles[i].x<<", "<<particles[i].y<<") @ "<<particles[i].theta<<" w: "<<particles[i].weight<<endl;
  }
 
  
  
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
/*
void display(){
  for(int i = 0; i<num_particles; ++i){
    std::cout<<"Particle #"<<i;
    std::cout<<"\tx: "<<particles[i].x<<"\ty: "<<particles[i].y<<"\tq: "<<particles[i].theta<<"\tw: "<<particles[i].weight<<std::endl;
  }
}
*/