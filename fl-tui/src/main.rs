use fl_sim::Simulation;
use nalgebra::vector;

fn main() {
  let mut simulation = Simulation::new(
    vector![10.0, 10.0, 10.0],
    fl_sim::Settings {
      delta_time:       0.1,
      smoothing_length: 1.0,
      rest_density:     1.0,
      iterations:       10,
      constraint:       1.0,
      viscosity:        1.0,
    },
  );

  let mut first = true;

  loop {
    simulation.tick();

    if !first {
      print!("\x1B[22A");
    }

    for _ in 0..50 {
      print!("=");
    }
    println!();

    for _y in 0..20 {
      for _x in 0..50 {
        print!(".");
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
