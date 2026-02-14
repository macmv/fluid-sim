use nalgebra::{Point2, Vector2, point, vector};
use std::f32::consts::PI;

use crate::space::SpatialIndex;

mod space;

pub struct Settings {
  pub delta_time:       f32,
  pub smoothing_length: f32,
  pub rest_density:     f32,
  pub iterations:       u32,
  pub constraint:       f32,
  pub viscosity:        f32,
}

pub struct Simulation {
  settings:  Settings,
  size:      Vector2<f32>,
  particles: Vec<Particle>,
  index:     SpatialIndex,
}

#[derive(Debug)]
pub struct Particle {
  pub position: Point2<f32>,
  pub velocity: Vector2<f32>,
  pub density:  f32,

  density_lambda: f32,
  predicted:      Point2<f32>,
}

const GRAVITY: Vector2<f32> = vector![0.0, -9.8];
const REST_DENSITY: f32 = 1000.0; // kg/m^2
const PARTICLE_SPACING: f32 = 0.5; // particles/m
const PARTICLE_MASS: f32 = 250.0; // kg
const LAMBDA_EPSILON: f32 = 1e-6;
const ITERATIONS: u32 = 15;
const SCORR_K: f32 = 0.001;
const SCORR_N: i32 = 4;
const SCORR_Q: f32 = 0.3;

impl Simulation {
  pub fn new(size: Vector2<f32>, settings: Settings) -> Simulation {
    Simulation {
      settings,
      size,
      particles: vec![],
      index: SpatialIndex::new(size, 2.0 * PARTICLE_SPACING),
    }
  }

  pub fn add_particle(&mut self, pos: Point2<f32>) {
    self.particles.push(Particle {
      position: pos,
      velocity: vector![0.0, 0.0],
      density:  0.0,

      density_lambda: 0.0,
      predicted:      point![0.0, 0.0],
    });
  }

  pub fn apply_repulsion(&mut self, center: Point2<f32>, radius: f32, strength: f32) {
    for particle in self.particles.iter_mut() {
      let delta = particle.position - center;
      let distance = delta.norm();
      if distance == 0.0 || distance >= radius {
        continue;
      }

      let direction = delta / distance;
      let falloff = 1.0 - distance / radius;
      particle.velocity += direction * (strength * falloff * self.settings.delta_time);
    }
  }

  pub fn tick(&mut self) {
    for (id, particle) in self.particles.iter_mut().enumerate() {
      particle.velocity += GRAVITY * self.settings.delta_time;
      particle.predicted = particle.position + particle.velocity * self.settings.delta_time;

      if particle.predicted.x < 0.0 {
        particle.predicted.x = 0.0;
        particle.velocity.x *= -self.settings.constraint;
      } else if particle.predicted.x > self.size.x {
        particle.predicted.x = self.size.x;
        particle.velocity.x *= -self.settings.constraint;
      }

      if particle.predicted.y < 0.0 {
        particle.predicted.y = 0.0;
        particle.velocity.y *= -self.settings.constraint;
      } else if particle.predicted.y > self.size.y {
        particle.predicted.y = self.size.y;
        particle.velocity.y *= -self.settings.constraint;
      }

      self.index.move_particle(id as u32, particle.predicted);
    }

    for _ in 0..ITERATIONS {
      for id in 0..self.particles.len() {
        let mut estimated_density = 0.0;
        let mut gradient_sum = vector![0.0, 0.0];
        let mut gradient_sum_squared = 0.0;

        for neighbor in self.index.neighbors(id as u32) {
          let p = &self.particles[id];
          let n = &self.particles[neighbor as usize];
          let delta = p.predicted - n.predicted;
          let distance = delta.norm();

          if distance >= self.index.radius() {
            continue;
          }

          let weight = kernel_poly6(distance, self.index.radius());
          estimated_density += PARTICLE_MASS * weight;

          let gradient =
            kernel_spiky_gradient(delta, self.index.radius()) * (PARTICLE_MASS / REST_DENSITY);
          gradient_sum += gradient;
          gradient_sum_squared += gradient.norm_squared();
        }

        gradient_sum_squared += gradient_sum.norm_squared();
        let density_constraint = (estimated_density / REST_DENSITY - 1.0).max(-0.005);

        self.particles[id].density = estimated_density;
        self.particles[id].density_lambda =
          -density_constraint / (gradient_sum_squared + LAMBDA_EPSILON);
      }

      let mut position_deltas = vec![vector![0.0, 0.0]; self.particles.len()];
      for id in 0..self.particles.len() {
        let mut total_position_delta = vector![0.0, 0.0];

        for neighbor in self.index.neighbors(id as u32) {
          let p = &self.particles[id];
          let n = &self.particles[neighbor as usize];
          let delta = p.predicted - n.predicted;
          let distance = delta.norm();

          if distance >= self.index.radius() {
            continue;
          }

          let kernel_gradient = kernel_spiky_gradient(delta, self.index.radius());
          let s_corr = tensile_correction(distance, self.index.radius());

          // This makes the correction symmetric:
          // equal and opposite influence between i and j (momentum-friendly).
          total_position_delta +=
            (p.density_lambda + n.density_lambda + s_corr) * PARTICLE_MASS * kernel_gradient;
        }

        position_deltas[id] = total_position_delta / REST_DENSITY;
      }

      for id in 0..self.particles.len() {
        self.particles[id].predicted += position_deltas[id];
        self.particles[id].predicted.x = self.particles[id].predicted.x.clamp(0.0, self.size.x);
        self.particles[id].predicted.y = self.particles[id].predicted.y.clamp(0.0, self.size.y);
        self.index.move_particle(id as u32, self.particles[id].predicted);
      }
    }

    for particle in self.particles.iter_mut() {
      particle.velocity = (particle.predicted - particle.position) / self.settings.delta_time;
      particle.position = particle.predicted;
    }
  }

  pub fn particles(&self) -> impl Iterator<Item = &Particle> { self.particles.iter() }
}

fn kernel_poly6(distance: f32, radius: f32) -> f32 {
  const NUMERATOR_2D: f32 = 4.0;
  const FACTOR_2D: f32 = 1.0;

  if distance >= radius {
    return 0.0;
  }

  let coeff = NUMERATOR_2D / (FACTOR_2D * PI * radius.powi(8));
  coeff * (radius.powi(2) - distance.powi(2)).powi(3)
}

fn kernel_spiky_gradient(displacement: Vector2<f32>, radius: f32) -> Vector2<f32> {
  const NUMERATOR_2D: f32 = -30.0;
  const FACTOR_2D: f32 = 1.0;

  let distance = displacement.norm();
  if distance == 0.0 || distance >= radius {
    return vector![0.0, 0.0];
  }

  let coeff = NUMERATOR_2D / (FACTOR_2D * PI * radius.powi(5));
  (displacement / distance) * (coeff * (radius - distance).powi(2))
}

fn tensile_correction(distance: f32, radius: f32) -> f32 {
  let numerator = kernel_poly6(distance, radius);
  let denominator = kernel_poly6(SCORR_Q * radius, radius);
  if denominator <= 0.0 {
    return 0.0;
  }
  -SCORR_K * (numerator / denominator).powi(SCORR_N)
}
