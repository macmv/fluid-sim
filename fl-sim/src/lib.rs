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
}

impl Simulation {
  pub fn new(size: Vector3<f32>, settings: Settings) -> Simulation {
    Simulation { settings, size, particles: vec![], index: SpatialIndex::new(size, 1.0) }
  }

  pub fn add_particle(&mut self, pos: Point3<f32>) {
    self.particles.push(Particle {
      position: pos,
      velocity: vector![0.0, 0.0, 0.0],
      density:  0.0,
    });
  }

  pub fn tick(&mut self) {
    const GRAVITY: Vector3<f32> = vector![0.0, 9.8, 0.0];

    for (id, particle) in self.particles.iter_mut().enumerate() {
      particle.velocity += GRAVITY * self.settings.delta_time;
      particle.position += particle.velocity * self.settings.delta_time;

      if particle.position.y >= self.size.y {
        particle.position.y = self.size.y;
        particle.velocity.y *= -self.settings.constraint;
      }

      self.index.move_particle(id as u32, particle.position);
    }
  }

  pub fn particle_positions(&self) -> impl Iterator<Item = Point3<f32>> {
    self.particles.iter().map(|p| p.position)
  }
}

fn kernel_poly6(distance: f32, radius: f32) -> f32 {
  if distance >= radius {
    return 0.0;
  }

  // TODO: Normalization constant?
  (radius.powi(2) - distance.powi(2)).powi(3)
}

fn kernel_spiky_gradient(displacement: Vector3<f32>, radius: f32) -> Vector3<f32> {
  let distance = displacement.norm();
  if distance == 0.0 || distance >= radius {
    return vector![0.0, 0.0, 0.0];
  }

  (displacement / distance) * -(radius - distance).powi(2)
}
