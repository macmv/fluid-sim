use nalgebra::{Point3, Vector3, vector};

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
  size:      Vector3<f32>,
  particles: Vec<Particle>,
  index:     SpatialIndex,
}

struct Particle {
  position: Point3<f32>,
  velocity: Vector3<f32>,
  density:  f32,
  pressure: f32,
}

impl Simulation {
  pub fn new(size: Vector3<f32>, settings: Settings) -> Simulation {
    Simulation {
      settings,
      size,
      particles: vec![],
      index: SpatialIndex::new(size, vector![10, 10, 10]),
    }
  }

  pub fn tick(&mut self) {}

  pub fn particle_positions(&self) -> impl Iterator<Item = Point3<f32>> {
    self.particles.iter().map(|p| p.position)
  }
}
