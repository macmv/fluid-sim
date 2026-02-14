use eframe::egui;
use egui_plot::{Plot, PlotBounds, PlotPoints, Points};
use fl_sim::{Settings, Simulation};
use nalgebra::{point, vector};

const WORLD_WIDTH: f32 = 40.0;
const WORLD_HEIGHT: f32 = 20.0;
const MOUSE_FORCE_RADIUS: f32 = 2.0;
const MOUSE_FORCE_STRENGTH: f32 = 500.0;

struct App {
  simulation: Simulation,
  paused:     bool,
}

fn make_simulation() -> Simulation {
  let mut simulation = Simulation::new(
    vector![WORLD_WIDTH, WORLD_HEIGHT],
    Settings {
      delta_time:       0.01,
      smoothing_length: 1.0,
      rest_density:     1000.0,
      iterations:       10,
      constraint:       0.0,
      viscosity:        0.0,
    },
  );

  for y in 5..35 {
    for x in 20..50 {
      simulation.add_particle(point![x as f32 / 2.0, y as f32 / 2.0]);
    }
  }

  simulation
}

impl App {
  fn new() -> Self { Self { simulation: make_simulation(), paused: false } }
}

impl eframe::App for App {
  fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
    let mut single_step = false;
    ctx.input(|input| {
      if input.key_pressed(egui::Key::Space) {
        self.paused = !self.paused;
      }
      if input.key_pressed(egui::Key::S) {
        single_step = true;
      }
      if input.key_pressed(egui::Key::R) {
        self.simulation = make_simulation();
      }
    });

    let mut do_tick = !self.paused || single_step;

    egui::TopBottomPanel::top("controls").show(ctx, |ui| {
      ui.horizontal(|ui| {
        let state = if self.paused { "Paused" } else { "Running" };
        ui.label(format!("State: {state}"));
        if ui.button("Play/Pause [space]").clicked() {
          self.paused = !self.paused;
        }
        if ui.button("Single [s]").clicked() {
          do_tick = true;
        }
        if ui.button("Restart [r]").clicked() {
          self.simulation = make_simulation();
        }
      });
    });

    let pointer_down = ctx.input(|input| input.pointer.primary_down());
    let mut repel_center: Option<(f32, f32)> = None;

    egui::CentralPanel::default().show(ctx, |ui| {
      let particles = PlotPoints::from_iter(
        self.simulation.particles().map(|p| [p.position.x as f64, p.position.y as f64]),
      );

      let points = Points::new("particles", particles)
        .radius(2.0)
        .color(egui::Color32::from_rgb(80, 180, 255));

      Plot::new("fluid_particles")
        .allow_boxed_zoom(false)
        .allow_drag(false)
        .allow_scroll(false)
        .allow_zoom(false)
        .show_axes([true, true])
        .data_aspect(1.0)
        .show(ui, |plot_ui| {
          plot_ui.set_plot_bounds(PlotBounds::from_min_max(
            [0.0, 0.0],
            [WORLD_WIDTH as f64, WORLD_HEIGHT as f64],
          ));
          plot_ui.points(points);
          if pointer_down {
            if let Some(pointer) = plot_ui.pointer_coordinate() {
              repel_center = Some((pointer.x as f32, pointer.y as f32));
            }
          }
        });
    });

    if let Some((x, y)) = repel_center {
      self.simulation.apply_repulsion(point![x, y], MOUSE_FORCE_RADIUS, MOUSE_FORCE_STRENGTH);
      do_tick = true;
    }

    if do_tick {
      self.simulation.tick();
    }

    ctx.request_repaint();
  }
}

fn main() -> eframe::Result<()> {
  let options = eframe::NativeOptions::default();
  eframe::run_native("fl-gui", options, Box::new(|_cc| Ok(Box::new(App::new()))))
}
