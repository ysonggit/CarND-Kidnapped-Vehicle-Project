/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "particle_filter.h"
#include <limits>
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 105;  // TODO: Set the number of particles

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for(int i=0; i < num_particles; i++){
        Particle p;
        p.id = i;
        p.x =  dist_x(gen);
        p.y =  dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1.0;
        particles.push_back(p);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    // gaussian noise generator
    default_random_engine gen;
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);
    for(int i = 0; i< particles.size(); i++){
        Particle p_0 = particles[i];
        if (fabs(yaw_rate) > 0.00001 ){
            // equations are given in lesson 6, section 8
            particles[i].x = p_0.x + (velocity / yaw_rate) * (sin(p_0.theta + yaw_rate * delta_t) - sin(p_0.theta));
            particles[i].y = p_0.y + (velocity / yaw_rate) * (cos(p_0.theta) - cos(p_0.theta + yaw_rate * delta_t));
            particles[i].theta = p_0.theta + yaw_rate * delta_t;
        }else{
            // move along the current direction
            particles[i].x = p_0.x + velocity * delta_t * cos(p_0.theta);
            particles[i].y = p_0.y + velocity * delta_t * sin(p_0.theta);
        }
        // add noise to particles position
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
    //cout<<"------------ prediction step ---------------- "<<endl;
    //printParticles();
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (LandmarkObs & cur_obs : observations){
        double min_dist = std::numeric_limits<double>::max();
        for(int i=0; i< predicted.size(); i++){
            double cur_dist = dist(cur_obs.x, cur_obs.y, predicted[i].x, predicted[i].y);
            if (cur_dist <= min_dist){
                min_dist = cur_dist;
                cur_obs.id = predicted[i].id; // predicted[i].id is NOT i!!!!!
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    //  Ah-ha! I knew this book well! Refer to the python solutionof quiz 16 in Lesson 6
    for(int i = 0; i< particles.size(); i++){
        Particle p = particles[i];
        vector<LandmarkObs> predictions_in_range;
        for(auto lm : map_landmarks.landmark_list){
            double lm_dist = dist(lm.x_f, lm.y_f, p.x, p.y);
            if (lm_dist <= sensor_range){
                LandmarkObs lm_obs;
                lm_obs.id = lm.id_i;
                lm_obs.x = lm.x_f;
                lm_obs.y = lm.y_f;
                predictions_in_range.push_back(lm_obs);
            }
        }
        // translate observations from vehicle coordinate to map coordinate
        vector<LandmarkObs> transformed_obs;
        for(auto obs : observations){
            double x_map = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);
            double y_map = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
            LandmarkObs trans_obs;
            trans_obs.id = obs.id;
            trans_obs.x = x_map;
            trans_obs.y = y_map;
            transformed_obs.push_back(trans_obs);
        }
        // compute weights
        dataAssociation(predictions_in_range, transformed_obs);
        particles[i].weight = 1.0;
        double coefficient = 1.0/(2 * M_PI * std_landmark[0] * std_landmark[1]);
        for(auto obs: transformed_obs){
            // x and y are the observations in map coordinates from landmarks
            double x = obs.x;
            double y = obs.y;
            // mu_x and mu_y are the coordinates of the nearest landmarks
            double mu_x = 0.0, mu_y = 0.0;
            for(auto pred: predictions_in_range){
                if (pred.id == obs.id){
                    mu_x = pred.x;
                    mu_y = pred.y;
                    break;
                }
            }
            double a = pow((x-mu_x),2)/(2*pow(std_landmark[0],2));
            double b = pow((y-mu_y),2)/(2*pow(std_landmark[1],2));
            // weight equal to multivariate Gaussia products
            particles[i].weight *= coefficient * exp(-(a+b));
            //cout<<"particle ["<<i<<"].weight = "<< particles[i].weight<<endl;
        }
    }
    //cout<<"------------ update weights step ---------------- "<<endl;
    //printParticles();
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    
    vector<Particle> sampled_particles;
    vector<double> weights;
    for(auto p : particles){
        weights.push_back(p.weight);
    }
    uniform_int_distribution<int> dist_int(0, num_particles-1); // uniformly distributed on the closed interval [a, b]

    auto index = dist_int(gen);
    double w_max = *max_element(weights.begin(), weights.end());
    uniform_real_distribution<double> dist_double (0.0, 2*w_max);
    double beta = 0.0;
    for (int i = 0; i<num_particles; i++){
        beta += dist_double(gen);
        while(weights[index] < beta){
            beta -= weights[index];
            index = (index+1) % num_particles;
        }
        sampled_particles.push_back(particles[index]);
    }
    particles = sampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

void ParticleFilter::printParticles(){
    for(auto p : particles){
        std::cout<< "Particle [" << p.id<<"]: ("<<p.x<<", "<<p.y<<", "<<p.theta<<"), weight: " << p.weight << std::endl;
    }
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
