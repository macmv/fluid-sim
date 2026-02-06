use std::collections::{HashMap, HashSet};

use nalgebra::{Vector3, vector};

pub struct SpatialIndex {
  cell_size:     Vector3<f32>,
  size:          Vector3<u32>,
  cells:         Vec<HashSet<ParticleId>>,
  reverse_cells: HashMap<ParticleId, u32>,
}

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
struct ParticleId(u32);

impl SpatialIndex {
  pub fn new(world_size: Vector3<f32>, size: Vector3<u32>) -> SpatialIndex {
    SpatialIndex {
      cell_size: vector![
        world_size.x / size.x as f32,
        world_size.y / size.y as f32,
        world_size.z / size.z as f32
      ],
      size,
      cells: vec![HashSet::new(); (size.x * size.y * size.z) as usize],
      reverse_cells: HashMap::new(),
    }
  }

  fn pos_to_cell(&self, pos: Vector3<f32>) -> u32 {
    if pos.x < 0.0 || pos.y < 0.0 || pos.z < 0.0 {
      panic!("position outside of index");
    }

    let x = (pos.x / self.cell_size.x) as u32;
    let y = (pos.y / self.cell_size.y) as u32;
    let z = (pos.z / self.cell_size.z) as u32;
    if x >= self.size.x || y >= self.size.y || z >= self.size.z {
      panic!("position outside of index");
    }

    let index = x + self.size.x * (y + self.size.y * z);
    index
  }

  pub fn move_particle(&mut self, id: u32, position: Vector3<f32>) {
    let new_cell = self.pos_to_cell(position);
    let prev_cell = self.reverse_cells.get(&ParticleId(id)).copied();
    if prev_cell != Some(new_cell) {
      self.reverse_cells.remove(&ParticleId(id));
      self.reverse_cells.insert(ParticleId(id), new_cell);

      if let Some(prev) = prev_cell {
        self.cells[prev as usize].remove(&ParticleId(id));
      }
      self.cells[new_cell as usize].insert(ParticleId(id));
    }
  }

  pub fn neighbors(&self, id: u32) -> impl Iterator<Item = u32> {
    let cell = self.reverse_cells.get(&ParticleId(id)).copied().expect("no such particle");
    let x = cell % self.size.x;
    let y = (cell / self.size.x) % self.size.y;
    let z = cell / (self.size.x * self.size.y);

    let x_range = x.saturating_sub(1)..=(x.saturating_add(1).min(self.size.x - 1));
    let y_range = y.saturating_sub(1)..=(y.saturating_add(1).min(self.size.y - 1));
    let z_range = z.saturating_sub(1)..=(z.saturating_add(1).min(self.size.z - 1));

    z_range.into_iter().flat_map(move |z| {
      let x_range = x_range.clone();
      y_range.clone().into_iter().flat_map(move |y| {
        x_range.clone().into_iter().flat_map(move |x| {
          let cell = x + self.size.x * (y + self.size.y * z);
          self.cells[cell as usize].iter().map(|p| p.0)
        })
      })
    })
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn it_works() {
    let mut index = SpatialIndex::new(vector![10.0, 10.0, 10.0], vector![10, 10, 10]);
    index.move_particle(0, vector![0.5, 0.5, 0.5]);

    assert_eq!(index.cells[0].len(), 1);

    index.move_particle(0, vector![1.5, 0.5, 0.5]);

    assert_eq!(index.cells[0].len(), 0);
    assert_eq!(index.cells[1].len(), 1);
  }
}
