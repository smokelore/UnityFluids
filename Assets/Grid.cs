using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cell {
	public int[] index = new int[3];
	public Vector3 position;
	public Vector3 velocity;
	public Vector3 prev_velocity;
	public float density;
	public float prev_density;
	public bool isBoundary = false;
	public ParticleSystem.Particle particle;

	public Vector3 allowVelocity = Vector3.one;

	public void SetVelocity(Vector3 new_velocity) {
		this.prev_velocity = this.velocity;
		new_velocity = new Vector3(allowVelocity.x * new_velocity.x, allowVelocity.y + new_velocity.y, allowVelocity.z + new_velocity.z);	// only  allow changes to allowed velocities (because boundaries)
		this.velocity = new_velocity;		
	}

	public void SetBoundary(Vector3 allowVelocity) {
		allowVelocity = new Vector3(Mathf.Clamp(Mathf.Round(allowVelocity.x), 0f, 1f), Mathf.Clamp(Mathf.Round(allowVelocity.y), 0f, 1f), Mathf.Clamp(Mathf.Round(allowVelocity.z), 0f, 1f));
		this.allowVelocity = allowVelocity;
		if (allowVelocity != Vector3.one) {
			this.isBoundary = true;
		} else {
			this.isBoundary = false;
		}
		//Debug.Log("" + index[0] + " " + index[1] + " " + index[2] + " : " + isBoundary);
	}

	public void SetDensity(float new_density) {
		new_density = Mathf.Clamp(new_density, 0f, 1f);
		this.prev_density = this.density;
		this.density = new_density;

		if (!isBoundary) {
			ParticleSystem.Particle p = particle;
			p.color = Color.Lerp(Grid.instance.lowColor, Grid.instance.highColor, new_density);
			particle = p;
		}
	}

	public void ChangeDensity(float delta_density) {
		SetDensity(this.density + delta_density);
	}

	public bool PointIsWithin(Vector3 testPosition) {
		if (!isBoundary) {
			Vector3 minExtents = this.position - Vector3.one * this.particle.size/2;
			Vector3 maxExtents = this.position + Vector3.one * this.particle.size/2;
			return (testPosition.x < maxExtents.x && testPosition.x >= minExtents.x
					&& testPosition.y < maxExtents.y && testPosition.y >= minExtents.y
					&& testPosition.z < maxExtents.z && testPosition.z >= minExtents.z);
		} else {
			return false;
		}
	}
}

public class Source {
	public int[] index = new int[3];
	public float amount;

	private float initAmount;
	private float initTime;
	private float duration = 2f;

	public Source(int[] index, float amount) {
		this.index = index;
		this.amount = amount;
		this.initAmount = amount;
		this.initTime = Time.time;
	}

	public void UpdateEmission() {
		this.amount = Mathf.Lerp(Time.time - initTime, duration, initAmount);
	}

	public bool isAlive() {
		return (Time.time - initTime < duration);
	}
}

public class Grid : MonoBehaviour {
	public static Grid instance;
	public ParticleSystem system;
	public Vector3 GRID_RESOLUTION; 
	public Cell[,,] cells;

	private Vector3 ORIGIN;
	public float GRID_SIZE;
	private int[] RESOLUTION = new int[3];
	private float PARTICLE_SIZE;
	public Color highColor, lowColor;

	private float dt = 0;
	public float DIFFUSION_RATE;

	public List<Source> sources = new List<Source>();
	private ParticleSystem.Particle[] particles;

	public GameObject UserBlock;

	public int relaxationIterations;

	void Start() {
		instance = this;
		this.system = this.GetComponent<ParticleSystem>();
		this.RESOLUTION = new int[3] {(int) (GRID_RESOLUTION.x+2), (int) (GRID_RESOLUTION.y+2), (int) (GRID_RESOLUTION.z+2)};
		this.cells = new Cell[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		this.particles = new ParticleSystem.Particle[(RESOLUTION[0]-2) * (RESOLUTION[1]-2) * (RESOLUTION[2]-2)];
		this.ORIGIN = this.transform.position;
		
		this.PARTICLE_SIZE = GRID_SIZE / Mathf.Max(GRID_RESOLUTION.x, GRID_RESOLUTION.y, GRID_RESOLUTION.z);
		InitializeGrid();
	}

	void Update() {
		this.dt = Time.deltaTime;
		DetectInput();
		Emit();
		Diffuse();
		UpdateParticleSystem();
	}

	public void InitializeGrid() {
		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					Vector3 pos = ORIGIN + new Vector3(i, j, k) * PARTICLE_SIZE;
					cells[i, j, k] = new Cell();
					cells[i, j, k].position = pos;
					cells[i, j, k].index = new int[] {i, j, k};

					Vector3 allowVelocity = Vector3.one;
					if (i == 0 || i == RESOLUTION[0]-1) {
						// boundary front, back
						allowVelocity.z = 0f;
					}

					if (j == 0 || j == RESOLUTION[1]-1) {
						// boundary top, bottom
						allowVelocity.y = 0f;
					} 

					if (k == 0 || k == RESOLUTION[2]-1) {	
						// boundary left, right
						allowVelocity.x = 0f;
					}

					cells[i, j, k].SetBoundary(allowVelocity);

					if (!cells[i, j, k].isBoundary) {
						ParticleSystem.Particle p = new ParticleSystem.Particle();
						p.size = PARTICLE_SIZE;
						p.position = pos;
						cells[i, j, k].particle = p;
						cells[i, j, k].SetDensity(Random.Range(0f, 1f)); //TEMPORARY, make just 0f later
						UpdateParticle(i, j, k);
					}
				}
			}
		}
	}

	public void UpdateParticleSystem() {
		system.SetParticles(particles, particles.Length);
	}

	public void UpdateParticle(int i, int j, int k) {
		particles[(i-1) + (j-1) * (RESOLUTION[0]-2) + (k-1) * (RESOLUTION[1]-2) * (RESOLUTION[2]-2)] = cells[i, j, k].particle;
	}

	/*public ParticleSystem.Particle[] GetFlatParticles() {
		ParticleSystem.Particle[] flattened = new ParticleSystem.Particle[(RESOLUTION[0]-2) * (RESOLUTION[1]-2) * (RESOLUTION[2]-2)];
		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					flattened[(i-1) + (j-1) * (RESOLUTION[0]-2) + (k-1) * (RESOLUTION[1]-2) * (RESOLUTION[2]-2)] = cells[i,j,k].particle;
				}
			}
		}
		
		//for (int i = 0; i < flattened.Length; i++) {
		//	Debug.Log("" + i + " : " + flattened[i].position);
		//}

		return flattened;
	}*/

	public ParticleSystem.Particle GetParticle(int i, int j, int k) {
		return cells[i, j, k].particle;
	}

	public void SetParticle(int i, int j, int k, ParticleSystem.Particle p) {
		cells[i, j, k].particle = p;
	}

	public void DetectInput() {
		int i = (int) Mathf.Round((UserBlock.transform.position.x - ORIGIN.x) / (PARTICLE_SIZE));
		int j = (int) Mathf.Round((UserBlock.transform.position.y - ORIGIN.y) / (PARTICLE_SIZE));
		int k = (int) Mathf.Round((UserBlock.transform.position.z - ORIGIN.z) / (PARTICLE_SIZE));
		//Debug.Log("" + i + " " + j + " " + k);
		if (i > 0 && i < RESOLUTION[0]-1 && j > 0 && j < RESOLUTION[1]-1 && k > 0 && k < RESOLUTION[2]-1) {
			sources.Add(new Source(new int[] {i, j, k}, 3f));
			//Debug.Log("within cell " + i + " " + j + " " + k);
		}
	}

	public void Emit() {
		for (int a = 0; a < sources.Count; a++) {
			if (!sources[a].isAlive()) {
				sources.RemoveAt(a);
			} else {
				sources[a].UpdateEmission();
			}
		}

		foreach (Source s in sources) {
			int i = s.index[0];
			int j = s.index[1];
			int k = s.index[2];
			cells[i, j, k].ChangeDensity(s.amount * dt);
			UpdateParticle(i, j, k);	
		}

		//	Debug.Log(sources.Count);
	}

	public void Diffuse() {
		float a = dt * DIFFUSION_RATE * RESOLUTION[0] * RESOLUTION[1] * RESOLUTION[2];

		float[,,] x0 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					x0[i, j, k] = cells[i, j, k].density;
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						float newDensity = (x0[i, j, k] + a * (cells[i+1, j, k].density + cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density + cells[i, j, k-1].density)) / (1+6*a);
						
						cells[i, j, k].SetDensity(newDensity); 
						UpdateParticle(i, j, k);
					}
				}
			}
		}
	}
}