use nalgebra::{Point3, Vector3};

pub struct Simulation {
  delta_time:       f32,
  smoothing_length: f32,
  rest_density:     f32,
  iterations:       u32,
  constraint:       f32,
  viscosity:        f32,

  size:      Vector3<f32>,
  particles: Vec<Particle>,
}

struct Particle {
  position: Point3<f32>,
  velocity: Vector3<f32>,
  density:  f32,
  pressure: f32,
}
