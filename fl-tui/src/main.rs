use fl_sim::Simulation;
use nalgebra::{point, vector};

fn main() {
  let mut simulation = Simulation::new(
    vector![50.0, 20.0, 1.0],
    fl_sim::Settings {
      delta_time:       0.1,
      smoothing_length: 1.0,
      rest_density:     1.0,
      iterations:       10,
      constraint:       0.5,
      viscosity:        1.0,
    },
  );

  for y in 0..10 {
    for x in 0..10 {
      simulation.add_particle(point![x as f32, y as f32, 0.0]);
    }
  }

  let mut first = true;
  let mut density = vec![vec![0.0; 50]; 20];

  loop {
    simulation.tick();
    for it in density.iter_mut().flatten() {
      *it = 0.0;
    }
    for pos in simulation.particle_positions() {
      density[pos.y.clamp(0.0, 19.0) as usize][pos.x.clamp(0.0, 49.0) as usize] += 0.5;
    }

    if !first {
      print!("\x1B[22A");
    }

    for _ in 0..50 {
      print!("=");
    }
    println!();

    for y in 0..20 {
      for x in 0..50 {
        let d = density[y][x];
        print!(
          "{}",
          match d {
            v if v >= 1.0 => '#',
            v if v >= 0.5 => '.',
            _ => ' ',
          }
        );
      }
      println!();
    }

    for _ in 0..50 {
      print!("=");
    }
    println!();

    first = false;

    std::thread::sleep(std::time::Duration::from_millis(100));
  }
}
