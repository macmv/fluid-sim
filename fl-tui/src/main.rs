use fl_sim::Simulation;
use nalgebra::{point, vector};

fn density_to_rgb(density: f32) -> (u8, u8, u8) {
  if density <= 1000.0 {
    let t = (density / 1000.0).clamp(0.0, 1.0);
    let r = 0.0;
    let g = 255.0 * t;
    let b = 255.0 * (1.0 - t);
    (r as u8, g as u8, b as u8)
  } else {
    let t = ((density - 1000.0) / 10_000.0).clamp(0.0, 1.0);
    let r = 255.0 * t;
    let g = 255.0 * (1.0 - t);
    let b = 0.0;
    (r as u8, g as u8, b as u8)
  }
}

fn main() {
  let mut simulation = Simulation::new(
    vector![50.0, 20.0],
    fl_sim::Settings {
      delta_time:       0.1,
      smoothing_length: 1.0,
      rest_density:     1.0,
      iterations:       10,
      constraint:       0.5,
      viscosity:        1.0,
    },
  );

  for y in 0..20 {
    for x in 20..30 {
      simulation.add_particle(point![x as f32 / 2.0, y as f32 / 2.0]);
    }
  }

  let mut first = true;
  let mut density = vec![vec![0.0; 50]; 20];

  loop {
    simulation.tick();
    for it in density.iter_mut().flatten() {
      *it = 0.0;
    }
    for p in simulation.particles() {
      density[p.position.y.clamp(0.0, 19.0) as usize][p.position.x.clamp(0.0, 49.0) as usize] +=
        p.density;
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
        let (r, g, b) = density_to_rgb(d);
        if d >= 0.1 {
          print!(
            "\x1b[38;2;{r};{g};{b}m{}",
            match () {
              _ if d >= 2.0 => '#',
              _ => '.',
            }
          );
        } else {
          print!(" ");
        }
      }
      println!("\x1b[0m");
    }

    for _ in 0..50 {
      print!("=");
    }
    println!();

    first = false;

    std::thread::sleep(std::time::Duration::from_millis(100));
  }
}
