use std::collections::HashSet;

use nalgebra::{Point2, Vector2, vector};

pub struct SpatialIndex<const N: usize> {
  radius:        f32,
  size:          Vector2<u32>,
  cells:         Vec<HashSet<ParticleId>>,
  reverse_cells: [u32; N], // index is particle id
}

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
struct ParticleId(u32);

impl<const N: usize> SpatialIndex<N> {
  pub fn new(world_size: Vector2<f32>, radius: f32) -> Self {
    let size =
      vector![(world_size.x / radius).ceil() as u32, (world_size.y / radius).ceil() as u32];

    SpatialIndex {
      radius,
      size,
      cells: vec![HashSet::new(); (size.x * size.y) as usize],
      reverse_cells: [0; N],
    }
  }

  pub fn radius(&self) -> f32 { self.radius }

  fn pos_to_cell(&self, pos: Point2<f32>) -> Option<u32> {
    if pos.x < 0.0 || pos.y < 0.0 {
      return None;
    }

    let x = (pos.x / self.radius) as u32;
    let y = (pos.y / self.radius) as u32;
    if x >= self.size.x || y >= self.size.y {
      return None;
    }

    Some(x + self.size.x * y)
  }

  pub fn move_particle(&mut self, id: u32, position: Point2<f32>) {
    if let Some(new_cell) = self.pos_to_cell(position) {
      let prev_cell = self.reverse_cells.get(id as usize).copied();
      if prev_cell != Some(new_cell) {
        self.reverse_cells[id as usize] = new_cell;

        if let Some(prev) = prev_cell {
          self.cells[prev as usize].remove(&ParticleId(id));
        }
        self.cells[new_cell as usize].insert(ParticleId(id));
      }
    }
  }

  pub fn neighbors(&self, id: u32) -> impl Iterator<Item = u32> {
    let iter = self.reverse_cells.get(id as usize).copied().map(|cell| {
      let x = cell % self.size.x;
      let y = cell / self.size.x;

      let x_range = x.saturating_sub(1)..=(x.saturating_add(1).min(self.size.x - 1));
      let y_range = y.saturating_sub(1)..=(y.saturating_add(1).min(self.size.y - 1));

      y_range.into_iter().flat_map(move |y| {
        x_range.clone().into_iter().flat_map(move |x| {
          let cell = x + self.size.x * y;
          self.cells[cell as usize].iter().map(|p| p.0)
        })
      })
    });

    iter.into_iter().flatten()
  }
}

#[cfg(test)]
mod tests {
  use nalgebra::point;

  use super::*;

  #[test]
  fn it_works() {
    let mut index = SpatialIndex::<1>::new(vector![10.0, 10.0], 1.0);
    index.move_particle(0, point![0.5, 0.5]);

    assert_eq!(index.cells[0].len(), 1);

    index.move_particle(0, point![1.5, 0.5]);

    assert_eq!(index.cells[0].len(), 0);
    assert_eq!(index.cells[1].len(), 1);
  }
}
