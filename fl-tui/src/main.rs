use fl_sim::Simulation;
use nalgebra::{point, vector};

fn density_to_rgb(density: f32) -> (u8, u8, u8) {
  let r;
  let g;
  let b;

  if density <= 1.0 {
    let t = density.clamp(0.0, 1.0);
    r = 0.0;
    g = 255.0 * t;
    b = 255.0 * (1.0 - t);
  } else {
    let t = (density - 1.0).clamp(0.0, 1.0);
    r = 255.0 * t;
    g = 255.0 * (1.0 - t);
    b = 0.0;
  }
  (r as u8, g as u8, b as u8)
}

fn main() {
  const WIDTH: usize = 50;
  const HEIGHT: usize = 20;

  let mut simulation = Simulation::new(vector![WIDTH as f32, HEIGHT as f32]);

  for y in 20..40 {
    for x in 20..80 {
      simulation.add_particle(point![x as f32 / 2.0, y as f32 / 2.0]);
    }
  }

  let mut first = true;
  let mut density = vec![vec![0.0; WIDTH]; HEIGHT];
  let mut counts = vec![vec![0u32; WIDTH]; HEIGHT];

  let mut next = std::time::Instant::now() + std::time::Duration::from_millis(10);

  loop {
    simulation.tick();
    for it in density.iter_mut().flatten() {
      *it = 0.0;
    }
    for it in counts.iter_mut().flatten() {
      *it = 0;
    }
    for p in simulation.particles() {
      let y = p.position.y.clamp(0.0, HEIGHT as f32 - 1.0) as usize;
      let x = p.position.x.clamp(0.0, WIDTH as f32 - 1.0) as usize;
      density[y][x] += p.density;
      counts[y][x] += 1;
    }
    for y in 0..HEIGHT {
      for x in 0..WIDTH {
        if counts[y][x] > 0 {
          density[y][x] /= counts[y][x] as f32;
        }
      }
    }

    if !first {
      print!("\x1B[22A");
    }

    for _ in 0..WIDTH {
      print!("=");
    }
    println!();

    for y in (0..HEIGHT).rev() {
      for x in 0..WIDTH {
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

    for _ in 0..WIDTH {
      print!("=");
    }
    println!();

    first = false;

    let now = std::time::Instant::now();
    if let Some(sleep) = next.checked_duration_since(now) {
      std::thread::sleep(sleep);
      next += std::time::Duration::from_millis(10);
    } else {
      next = now;
    }
  }
}
