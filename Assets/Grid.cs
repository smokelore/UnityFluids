using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cell {
	public int[] index = new int[3];
	public Vector3 position;
	//public Vector3 force;
	public Vector3 velocity;
	public Vector3 prev_velocity;
	public float density;
	public float prev_density;
	public bool isBoundary = false;
	public ParticleSystem.Particle particle;

	public Vector3 allowVelocity = Vector3.one;

	public void SetVelocity(Vector3 new_velocity) {
		this.prev_velocity = this.velocity;
		new_velocity = Vector3.ClampMagnitude(new_velocity, 1f);
		this.velocity = new Vector3(allowVelocity.x * new_velocity.x, allowVelocity.y * new_velocity.y, allowVelocity.z * new_velocity.z);	// only  allow changes to allowed velocities (because boundaries)
		//Debug.DrawLine(position, position + velocity, Color.red);
	}

	public void ForceVelocity(Vector3 new_velocity) {
		this.prev_velocity = this.velocity;
		new_velocity = Vector3.ClampMagnitude(new_velocity, 1f);
		//Debug.DrawLine(position, position + velocity, Color.red);
	}

	public void ChangeVelocity(Vector3 delta) {
		SetVelocity(1f * this.velocity + 1f * delta);
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
	public Vector3 direction;
	public float amount;

	public float initAmount;
	private float initTime;
	private float duration = 0.1f;

	public Source(int[] index, float amount) {
		Initialize(index, Vector3.zero, amount);
	}

	public Source(int[] index, Vector3 direction, float amount) {
		Initialize(index, direction, amount);
	}

	public Source(int i, int j, int k, float amount) {
		Initialize(new int[3] {i, j, k}, Vector3.zero, amount);
	}

	public Source(int i, int j, int k, Vector3 direction, float amount) {
		Initialize(new int[3] {i, j, k}, direction, amount);
	}

	private void Initialize(int[] index, Vector3 direction, float amount) {
		this.index = index;
		this.direction = Vector3.ClampMagnitude(direction, 1f);
		this.initAmount = amount;
		this.initTime = Time.time;

		Grid.instance.cells[index[0], index[1], index[2]].ChangeDensity(amount * Grid.instance.dt);
		Grid.instance.cells[index[0], index[1], index[2]].ChangeVelocity(amount * direction * Grid.instance.dt);
		//Grid.instance.cells[index[0], index[1], index[2]].ChangeVelocity(amount * direction);
	}

	public void UpdateEmission() {
		float progress = (Time.time - initTime)/duration;
		this.amount = Mathf.Lerp(0f, 1f, Mathf.Pow(progress+0.25f, 2f));
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

	public float dt = 0;
	public float DIFFUSION_RATE;

	public List<Source> sources = new List<Source>();
	private ParticleSystem.Particle[] particles;

	public GameObject UserObject;

	public int relaxationIterations;

	void Start() {
		instance = this;
		this.system = this.GetComponent<ParticleSystem>();
		this.RESOLUTION = new int[3] {(int) (GRID_RESOLUTION.x+2), (int) (GRID_RESOLUTION.y+2), (int) (GRID_RESOLUTION.z+2)};
		this.cells = new Cell[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		this.particles = new ParticleSystem.Particle[(RESOLUTION[0]) * (RESOLUTION[1]) * (RESOLUTION[2])];
		this.ORIGIN = this.transform.position;
		
		this.PARTICLE_SIZE = GRID_SIZE / Mathf.Max(GRID_RESOLUTION.x, GRID_RESOLUTION.y, GRID_RESOLUTION.z);
		InitializeGrid();
	}

	void Update() {
		this.dt = Time.deltaTime;

		DetectInput();

		Emit();

		Diffuse();

		Advect();

		Project();

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
					cells[i, j, k].SetVelocity(new Vector3(0, 0, 0));
					//cells[i, j, k].SetVelocity(new Vector3(Random.Range(-1f, 1f), Random.Range(-1f, 1f), Random.Range(-1f, 1f))/10f);

					//if (!cells[i, j, k].isBoundary) {
						ParticleSystem.Particle p = new ParticleSystem.Particle();
						p.size = PARTICLE_SIZE;
						p.position = pos;
						cells[i, j, k].particle = p;
						cells[i, j, k].SetDensity(0f);
						UpdateParticle(i, j, k);
					//}
				}
			}
		}
	}

	public void UpdateParticleSystem() {
		system.SetParticles(particles, particles.Length);
	}

	public void UpdateParticle(int i, int j, int k) {
		particles[(i) + (j) * (RESOLUTION[0]) + (k) * (RESOLUTION[1]) * (RESOLUTION[2])] = cells[i, j, k].particle;
		if (cells[i, j, k].isBoundary) {
			Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
		} else {
			Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.red);
		}
	}

	public ParticleSystem.Particle GetParticle(int i, int j, int k) {
		return cells[i, j, k].particle;
	}

	public void SetParticle(int i, int j, int k, ParticleSystem.Particle p) {
		cells[i, j, k].particle = p;
	}

	public void DetectInput() {
		for (int x = 0; x < RESOLUTION[0]; x++) {
			for (int y = 0; y < RESOLUTION[1]; y++) {
				for (int z = 0; z < RESOLUTION[2]; z++) {
					cells[x, y, z].prev_velocity = Vector3.zero;
					cells[x, y, z].prev_density = 0f;
				}
			}
		}

		int i = Mathf.RoundToInt((UserObject.transform.position.x - ORIGIN.x) / (PARTICLE_SIZE));
		int j = Mathf.RoundToInt((UserObject.transform.position.y - ORIGIN.y) / (PARTICLE_SIZE));
		int k = Mathf.RoundToInt((UserObject.transform.position.z - ORIGIN.z) / (PARTICLE_SIZE));
		//Debug.Log("" + i + " " + j + " " + k);
		if (i > 0 && i < RESOLUTION[0]-1 && j > 0 && j < RESOLUTION[1]-1 && k > 0 && k < RESOLUTION[2]-1) {
			Vector3 userVelocity = UserObject.GetComponent<CharacterController>().velocity;
			sources.Add(new Source(new int[] {i, j, k}, userVelocity, 10f * userVelocity.magnitude));
			//Debug.Log("within cell " + i + " " + j + " " + k);
			//Debug.Log(userVelocity);
		}
	}

	public void Emit() {
		foreach (Source s in sources) {
			int i = s.index[0];
			int j = s.index[1];
			int k = s.index[2];
			s.UpdateEmission();
			UpdateParticle(i, j, k);	
		}

		for (int a = 0; a < sources.Count; a++) {
			if (!sources[a].isAlive()) {
				sources.RemoveAt(a);
			} else {
				sources[a].UpdateEmission();
			}
		}
		//	Debug.Log(sources.Count);
	}

	public void Diffuse() {
		float a = dt * DIFFUSION_RATE * RESOLUTION[0] * RESOLUTION[1] * RESOLUTION[2];

		float[,,] d0 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		Vector3[,,] v0 = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					d0[i, j, k] = cells[i, j, k].density;
					v0[i, j, k] = cells[i, j, k].velocity;
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						float newDensity = (d0[i, j, k] + a * (cells[i+1, j, k].density + cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density + cells[i, j, k-1].density)) / (1+6*a);
						cells[i, j, k].SetDensity(newDensity);

						Vector3 newVelocity = (v0[i, j, k] + a * (cells[i+1, j, k].velocity + cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity + cells[i, j, k-1].velocity)) / (1+6*a);
						cells[i, j, k].SetVelocity(newVelocity); 
						
						UpdateParticle(i, j, k);
					}
				}
			}
			SetBoundaries();
		}
	}

	public void Advect() {
		float dt0i = dt * RESOLUTION[0];
		float dt0j = dt * RESOLUTION[1];
		float dt0k = dt * RESOLUTION[2];

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					float x = Mathf.Clamp(i - dt0i * cells[i, j, k].velocity.x, 0.5f, RESOLUTION[0]-1.5f);
					float y = Mathf.Clamp(j - dt0j * cells[i, j, k].velocity.y, 0.5f, RESOLUTION[1]-1.5f);
					float z = Mathf.Clamp(k - dt0k * cells[i, j, k].velocity.z, 0.5f, RESOLUTION[2]-1.5f);

					int i0 = (int) (x); 	int i1 = i0+1;
					int j0 = (int) (y); 	int j1 = j0+1;
					int k0 = (int) (z); 	int k1 = k0+1;

					float s1 = x - i0;	float s0 = 1-s1;
					float t1 = y - j0;	float t0 = 1-t1;
					float u1 = z - k0;	float u0 = 1-u1;

					//Debug.Log(i0 + " " + j0 + " " + k0);
					float d000 = cells[i0, j0, k0].prev_density;
					float d100 = cells[i1, j0, k0].prev_density;
					float d010 = cells[i0, j1, k0].prev_density;
					float d110 = cells[i1, j1, k0].prev_density;
					float d001 = cells[i0, j0, k1].prev_density;
					float d101 = cells[i1, j0, k1].prev_density;
					float d011 = cells[i0, j1, k1].prev_density;
					float d111 = cells[i1, j1, k1].prev_density;

					Vector3 v000 = cells[i0, j0, k0].prev_velocity;
					Vector3 v100 = cells[i1, j0, k0].prev_velocity;
					Vector3 v010 = cells[i0, j1, k0].prev_velocity;
					Vector3 v110 = cells[i1, j1, k0].prev_velocity;
					Vector3 v001 = cells[i0, j0, k1].prev_velocity;
					Vector3 v101 = cells[i1, j0, k1].prev_velocity;
					Vector3 v011 = cells[i0, j1, k1].prev_velocity;
					Vector3 v111 = cells[i1, j1, k1].prev_velocity;

					float newDensity = s0 * (t0 * (u0 * d000 + u1 * d001) + t1 * (u0 * d010 + u1 * d011)) + s1 * (t0 * (u0 * d100 + u1 * d101) + t1 * (u0 * d110 + u1 * d111));
				 	cells[i, j, k].SetDensity(newDensity);

				 	Vector3 newVelocity = s0 * (t0 * (u0 * v000 + u1 * v001) + t1 * (u0 * v010 + u1 * v011)) + s1 * (t0 * (u0 * v100 + u1 * v101) + t1 * (u0 * v110 + u1 * v111));
					cells[i, j, k].SetVelocity(newVelocity);

					UpdateParticle(i, j, k);
				}
			}
		}
	}

	public void Project() {
		Vector3 h = new Vector3(1f/RESOLUTION[0], 1f/RESOLUTION[1], 1f/RESOLUTION[2]);
		Vector3[,,] f0 = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		Vector3[,,] p = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					f0[i, j, k] = -0.5f * (h.x * (cells[i+1, j, k].velocity - cells[i-1, j, k].velocity) + h.y * (cells[i, j+1, k].velocity - cells[i, j-1, k].velocity) + h.z * (cells[i, j, k+1].velocity - cells[i, j, k-1].velocity));
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						p[i, j, k] = (f0[i, j, k] + p[i-1, j, k] + p[i+1, j, k] + p[i, j-1, k] + p[i, j+1, k] + p[i, j, k-1] + p[i, j, k+1]) / 6f;
					}
				}
			}
			SetBoundaries();
		}

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					Vector3 newVelocity = cells[i, j, k].velocity;
					newVelocity.x -= 0.5f * (p[i+1, j, k].x - p[i-1, j, k].x) / h.x;
					newVelocity.y -= 0.5f * (p[i, j+1, k].y - p[i, j-1, k].y) / h.y;
					newVelocity.z -= 0.5f * (p[i, j, k+1].z - p[i, j, k-1].z) / h.z;

					cells[i, j, k].SetVelocity(newVelocity);

					UpdateParticle(i, j, k);
				}
			}
		}
	}

	public void SetBoundaries() {
		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					bool left = (i == 0);
					bool bottom = (j == 0);
					bool front = (k == 0);

					bool right = (i == RESOLUTION[0]-1);
					bool top = (j == RESOLUTION[1]-1);
					bool back = (k == RESOLUTION[2]-1);
					if (left || right || bottom || front || top || back) {
						Vector3 newVelocity;
						float newDensity;

						if (left && bottom && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && bottom && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && top && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && top && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && bottom && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && bottom && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && top && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density + cells[i, j, k+1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && top && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k-1].velocity)/3f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density + cells[i, j, k-1].density)/3f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && bottom) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j+1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && top) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j-1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && front) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left && back) {
							newVelocity = (cells[i+1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i+1, j, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && bottom) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j+1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j+1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && top) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j-1, k].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j-1, k].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && front) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k+1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j, k+1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (right && back) {
							newVelocity = (cells[i-1, j, k].velocity + cells[i, j, k-1].velocity)/2f;
							cells[i, j, k].SetVelocity(newVelocity);
							newDensity = (cells[i-1, j, k].density + cells[i, j, k-1].density)/2f;
							cells[i, j, k].SetDensity(newDensity);
						} else if (left) {
							newVelocity.x = -1 *cells[i+1, j, k].velocity.x;
							newVelocity.y = 	cells[i+1, j, k].velocity.y;
							newVelocity.z = 	cells[i+1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i+1, j, k].density);
						} else if (right) {
							newVelocity.x = -1 *cells[i-1, j, k].velocity.x;
							newVelocity.y = 	cells[i-1, j, k].velocity.y;
							newVelocity.z = 	cells[i-1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i-1, j, k].density);
						} else if (bottom) {
							newVelocity.x = 	cells[i, j+1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j+1, k].velocity.y;
							newVelocity.z = 	cells[i, j+1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j+1, k].density);
						} else if (top) {
							newVelocity.x = 	cells[i, j-1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j-1, k].velocity.y;
							newVelocity.z = 	cells[i, j-1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j-1, k].density);
						} else if (front) {
							newVelocity.x = 	cells[i, j, k+1].velocity.x;
							newVelocity.y = 	cells[i, j, k+1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k+1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j, k+1].density);
						} else if (back) {
							newVelocity.x = 	cells[i, j, k-1].velocity.x;
							newVelocity.y =		cells[i, j, k-1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k-1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							cells[i, j, k].SetDensity(cells[i, j, k-1].density);
						}
						UpdateParticle(i, j, k);
						Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
					}
				}
			}
		}
	}
}