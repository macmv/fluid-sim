use criterion::{Criterion, criterion_group, criterion_main};
use fl_sim::Simulation;
use nalgebra::{point, vector};

const WORLD_WIDTH: f32 = 8.0;
const WORLD_HEIGHT: f32 = 8.0;
const PARTICLES: usize = 128;

fn make_simulation() -> Simulation<PARTICLES> {
  let mut simulation = Simulation::new(vector![WORLD_WIDTH, WORLD_HEIGHT]);

  let mut i = 0;
  for y in 4..4 + 8 {
    for x in 4..4 + 16 {
      simulation.set_particle(i, point![x as f32 / 4.0, y as f32 / 4.0]);
      i += 1;
    }
  }

  // simulation.add_barrier(point![4.0, 1.0], point![5.0, 2.0]);

  simulation
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut sim = make_simulation();

  c.bench_function("simulation tick", |b| b.iter(|| sim.tick()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
