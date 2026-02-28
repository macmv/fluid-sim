use alloc::vec::Vec;
use nalgebra::{Point2, Vector2, vector};

pub struct SpatialIndex<const N: usize> {
  radius:        f32,
  size:          Vector2<u8>,
  cells:         Vec<Set<ParticleId>>,
  reverse_cells: [u8; N], // index is particle id
}

#[derive(Clone)]
struct Set<T> {
  #[cfg(feature = "std")]
  contents: std::collections::HashSet<T>,
  #[cfg(not(feature = "std"))]
  contents: Vec<T>,
}

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
struct ParticleId(u8);

impl<const N: usize> SpatialIndex<N> {
  pub fn new(world_size: Vector2<f32>, radius: f32) -> Self {
    let size =
      vector![libm::ceilf(world_size.x / radius) as u8, libm::ceilf(world_size.y / radius) as u8];

    if size.x as u32 * size.y as u32 > u8::MAX as u32 + 1 {
      panic!("size is too large: {size:?}");
    }

    SpatialIndex {
      radius,
      size,
      cells: alloc::vec![Set::new(); size.x as usize * size.y as usize],
      reverse_cells: [0; N],
    }
  }

  pub fn radius(&self) -> f32 { self.radius }
  pub fn width(&self) -> u8 { self.size.x }
  pub fn height(&self) -> u8 { self.size.y }
  pub fn cell_count(&self) -> usize { self.cells.len() }

  fn pos_to_cell(&self, pos: Point2<f32>) -> Option<u8> {
    if pos.x < 0.0 || pos.y < 0.0 {
      return None;
    }

    let x = (pos.x / self.radius) as u8;
    let y = (pos.y / self.radius) as u8;
    if x >= self.size.x || y >= self.size.y {
      return None;
    }

    Some(x + self.size.x * y)
  }

  pub fn move_particle(&mut self, id: u8, position: Point2<f32>) {
    if let Some(new_cell) = self.pos_to_cell(position) {
      let prev_cell = self.reverse_cells.get(id as usize).copied();
      if prev_cell != Some(new_cell) {
        self.reverse_cells[id as usize] = new_cell;

        if let Some(prev) = prev_cell {
          self.cells[prev as usize].remove(ParticleId(id));
        }
        self.cells[new_cell as usize].insert(ParticleId(id));
      }
    }
  }

  pub fn neighbors(&self, id: u32) -> impl Iterator<Item = u8> {
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

impl<T> Set<T> {
  pub fn new() -> Self {
    Set {
      #[cfg(feature = "std")]
      contents:                              std::collections::HashSet::new(),
      #[cfg(not(feature = "std"))]
      contents:                              Vec::new(),
    }
  }

  pub fn iter(&self) -> impl Iterator<Item = &T> { self.contents.iter() }
}

impl<T: Eq + core::hash::Hash> Set<T> {
  pub fn insert(&mut self, value: T) {
    #[cfg(feature = "std")]
    self.contents.insert(value);
    #[cfg(not(feature = "std"))]
    self.contents.push(value);
  }

  pub fn remove(&mut self, value: T) {
    #[cfg(feature = "std")]
    self.contents.remove(&value);
    #[cfg(not(feature = "std"))]
    {
      if let Some(idx) = self.contents.iter().position(|v| *v == value) {
        self.contents.swap_remove(idx);
      }
    }
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

    assert_eq!(index.cells[0].contents.len(), 1);

    index.move_particle(0, point![1.5, 0.5]);

    assert_eq!(index.cells[0].contents.len(), 0);
    assert_eq!(index.cells[1].contents.len(), 1);
  }
}
