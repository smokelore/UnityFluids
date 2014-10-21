using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Cell {
	public int[] index = new int[3];
	public Vector3 position;
	//public Vector3 force;
	public Vector3 velocity;
	public Vector3 prev_velocity;
	public float[] density = new float[3];
	public float[] prev_density = new float[3];
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

	public void SetDensity(float new_density, int index) {
		new_density = Mathf.Clamp(new_density, 0f, 1f);
		this.prev_density[index] = this.density[index];
		this.density[index] = new_density;

		UpdateColor();
	}

	public void ChangeDensity(float delta_density, int index) {
		SetDensity(this.density[index] + delta_density, index);
	}

	public void UpdateColor() {
		Color newColor = new Color(density[0], density[1], density[2], Mathf.Max(density[0], density[1], density[2]));
		this.particle.color = newColor;
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
	public int colorIndex;

	public float initAmount;
	private float initTime;
	private float duration = 0.5f;

	public Source(int[] index, Vector3 direction, float amount, int colorIndex) {
		Initialize(index, direction, amount, colorIndex);
	}

	public Source(int i, int j, int k, float amount, int colorIndex) {
		Initialize(new int[3] {i, j, k}, Vector3.zero, amount, colorIndex);
	}

	public Source(int i, int j, int k, Vector3 direction, float amount, int colorIndex) {
		Initialize(new int[3] {i, j, k}, direction, amount, colorIndex);
	}

	private void Initialize(int[] index, Vector3 direction, float amount, int colorIndex) {
		this.index = index;
		this.direction = Vector3.ClampMagnitude(direction, 1f);
		this.initAmount = amount;
		this.initTime = Time.time;
		this.colorIndex = colorIndex;

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

	public float dt = 0;
	public float DIFFUSION_RATE;

	public List<Source> sources = new List<Source>();
	private ParticleSystem.Particle[] particles;

	public GameObject UserObject;

	public int relaxationIterations;

	public int chosenColor;

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
		SetBoundaries();

		Diffuse();
		SetBoundaries();

		Advect();
		SetBoundaries();

		Project();
		SetBoundaries();

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
						cells[i, j, k].SetDensity(0f, 0);
						cells[i, j, k].SetDensity(0f, 1);
						cells[i, j, k].SetDensity(0f, 2);
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
		int i = Mathf.RoundToInt((UserObject.transform.position.x - ORIGIN.x) / (PARTICLE_SIZE));
		int j = Mathf.RoundToInt((UserObject.transform.position.y - ORIGIN.y) / (PARTICLE_SIZE));
		int k = Mathf.RoundToInt((UserObject.transform.position.z - ORIGIN.z) / (PARTICLE_SIZE));
		//Debug.Log("" + i + " " + j + " " + k);
		if (i > 0 && i < RESOLUTION[0]-1 && j > 0 && j < RESOLUTION[1]-1 && k > 0 && k < RESOLUTION[2]-1) {
			Vector3 userVelocity = UserObject.GetComponent<CharacterController>().velocity;
			sources.Add(new Source(i, j, k, userVelocity, 100f * userVelocity.magnitude, chosenColor));
			//Debug.Log("within cell " + i + " " + j + " " + k);
			//Debug.Log(userVelocity);
		}

		if (Input.GetKeyUp("x")) {
			chosenColor++;
			while (chosenColor > 2) {
				chosenColor -= 3;
			}
		}
	}

	public void Emit() {
		foreach (Source s in sources) {
			int i = s.index[0];
			int j = s.index[1];
			int k = s.index[2];
			s.UpdateEmission();
			cells[i, j, k].ChangeDensity(s.amount * dt, s.colorIndex);
			cells[i, j, k].ChangeVelocity(s.amount * s.direction * dt);
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
		float[,,] d1 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		float[,,] d2 = new float[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];
		Vector3[,,] v0 = new Vector3[RESOLUTION[0], RESOLUTION[1], RESOLUTION[2]];

		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					d0[i, j, k] = cells[i, j, k].density[0];
					d1[i, j, k] = cells[i, j, k].density[1];
					d2[i, j, k] = cells[i, j, k].density[2];

					v0[i, j, k] = cells[i, j, k].velocity;
				}
			}
		}

		// Gauss-Seidel relaxation
		for (int r = 0; r < relaxationIterations; r++) {
			for (int i = 1; i < RESOLUTION[0]-1; i++) {
				for (int j = 1; j < RESOLUTION[1]-1; j++) {
					for (int k = 1; k < RESOLUTION[2]-1; k++) {
						float newDensity0 = (d0[i, j, k] + a * (cells[i+1, j, k].density[0] + cells[i-1, j, k].density[0] + cells[i, j+1, k].density[0] + cells[i, j-1, k].density[0] + cells[i, j, k+1].density[0] + cells[i, j, k-1].density[0])) / (1+6*a);
						cells[i, j, k].SetDensity(newDensity0, 0);

						float newDensity1 = (d1[i, j, k] + a * (cells[i+1, j, k].density[1] + cells[i-1, j, k].density[1] + cells[i, j+1, k].density[1] + cells[i, j-1, k].density[1] + cells[i, j, k+1].density[1] + cells[i, j, k-1].density[1])) / (1+6*a);
						cells[i, j, k].SetDensity(newDensity1, 1);

						float newDensity2 = (d2[i, j, k] + a * (cells[i+1, j, k].density[2] + cells[i-1, j, k].density[2] + cells[i, j+1, k].density[2] + cells[i, j-1, k].density[2] + cells[i, j, k+1].density[2] + cells[i, j, k-1].density[2])) / (1+6*a);
						cells[i, j, k].SetDensity(newDensity2, 2);

						Vector3 newVelocity = (v0[i, j, k] + a * (cells[i+1, j, k].velocity + cells[i-1, j, k].velocity + cells[i, j+1, k].velocity + cells[i, j-1, k].velocity + cells[i, j, k+1].velocity + cells[i, j, k-1].velocity)) / (1+6*a);
						cells[i, j, k].ChangeVelocity(newVelocity); 
						
						UpdateParticle(i, j, k);
					}
				}
			}
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
					float r000 = cells[i0, j0, k0].prev_density[0];
					float r100 = cells[i1, j0, k0].prev_density[0];
					float r010 = cells[i0, j1, k0].prev_density[0];
					float r110 = cells[i1, j1, k0].prev_density[0];
					float r001 = cells[i0, j0, k1].prev_density[0];
					float r101 = cells[i1, j0, k1].prev_density[0];
					float r011 = cells[i0, j1, k1].prev_density[0];
					float r111 = cells[i1, j1, k1].prev_density[0];

					float g000 = cells[i0, j0, k0].prev_density[1];
					float g100 = cells[i1, j0, k0].prev_density[1];
					float g010 = cells[i0, j1, k0].prev_density[1];
					float g110 = cells[i1, j1, k0].prev_density[1];
					float g001 = cells[i0, j0, k1].prev_density[1];
					float g101 = cells[i1, j0, k1].prev_density[1];
					float g011 = cells[i0, j1, k1].prev_density[1];
					float g111 = cells[i1, j1, k1].prev_density[1];

					float b000 = cells[i0, j0, k0].prev_density[2];
					float b100 = cells[i1, j0, k0].prev_density[2];
					float b010 = cells[i0, j1, k0].prev_density[2];
					float b110 = cells[i1, j1, k0].prev_density[2];
					float b001 = cells[i0, j0, k1].prev_density[2];
					float b101 = cells[i1, j0, k1].prev_density[2];
					float b011 = cells[i0, j1, k1].prev_density[2];
					float b111 = cells[i1, j1, k1].prev_density[2];

					Vector3 v000 = cells[i0, j0, k0].prev_velocity;
					Vector3 v100 = cells[i1, j0, k0].prev_velocity;
					Vector3 v010 = cells[i0, j1, k0].prev_velocity;
					Vector3 v110 = cells[i1, j1, k0].prev_velocity;
					Vector3 v001 = cells[i0, j0, k1].prev_velocity;
					Vector3 v101 = cells[i1, j0, k1].prev_velocity;
					Vector3 v011 = cells[i0, j1, k1].prev_velocity;
					Vector3 v111 = cells[i1, j1, k1].prev_velocity;

					float newDensity0 = s0 * (t0 * (u0 * r000 + u1 * r001) + t1 * (u0 * r010 + u1 * r011)) + s1 * (t0 * (u0 * r100 + u1 * r101) + t1 * (u0 * r110 + u1 * r111));
				 	cells[i, j, k].SetDensity(newDensity0, 0);

				 	float newDensity1 = s0 * (t0 * (u0 * g000 + u1 * g001) + t1 * (u0 * g010 + u1 * g011)) + s1 * (t0 * (u0 * g100 + u1 * g101) + t1 * (u0 * g110 + u1 * g111));
				 	cells[i, j, k].SetDensity(newDensity1, 1);

				 	float newDensity2 = s0 * (t0 * (u0 * b000 + u1 * b001) + t1 * (u0 * b010 + u1 * b011)) + s1 * (t0 * (u0 * b100 + u1 * b101) + t1 * (u0 * b110 + u1 * b111));
				 	cells[i, j, k].SetDensity(newDensity2, 2);

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
		}

		for (int i = 1; i < RESOLUTION[0]-1; i++) {
			for (int j = 1; j < RESOLUTION[1]-1; j++) {
				for (int k = 1; k < RESOLUTION[2]-1; k++) {
					Vector3 newVelocity = Vector3.zero;
					newVelocity.x = 0.5f * (p[i+1, j, k].x - p[i-1, j, k].x) / h.x;
					newVelocity.y = 0.5f * (p[i, j+1, k].y - p[i, j-1, k].y) / h.y;
					newVelocity.z = 0.5f * (p[i, j, k+1].z - p[i, j, k-1].z) / h.z;

					cells[i, j, k].SetVelocity(cells[i, j, k].velocity - newVelocity);

					UpdateParticle(i, j, k);
				}
			}
		}
	}

	public void SetBoundaries() {
		for (int i = 0; i < RESOLUTION[0]; i++) {
			for (int j = 0; j < RESOLUTION[1]; j++) {
				for (int k = 0; k < RESOLUTION[2]; k++) {
					if (cells[i, j, k].isBoundary) {
						bool left = (i == 0);
						bool bottom = (j == 0);
						bool front = (k == 0);

						bool right = (i == RESOLUTION[0]-1);
						bool top = (j == RESOLUTION[1]-1);
						bool back = (k == RESOLUTION[2]-1);

						Vector3 newVelocity;

						if ((left && top) || (left && bottom) || (right && top) || (right && bottom) || (left && front) || (left && back) || (right && front) || (right && back) || (top && front) || (top && back) || (bottom && front) || (bottom && back)) {
							// corners
							cells[i, j, k].SetVelocity(Vector3.zero);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (left) {
							newVelocity.x = -1 *cells[i+1, j, k].velocity.x;
							newVelocity.y = 	cells[i+1, j, k].velocity.y;
							newVelocity.z = 	cells[i+1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i+1, j, k].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (right) {
							newVelocity.x = -1 *cells[i-1, j, k].velocity.x;
							newVelocity.y = 	cells[i-1, j, k].velocity.y;
							newVelocity.z = 	cells[i-1, j, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i-1, j, k].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (bottom) {
							newVelocity.x = 	cells[i, j+1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j+1, k].velocity.y;
							newVelocity.z = 	cells[i, j+1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i, j+1, k].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (top) {
							newVelocity.x = 	cells[i, j-1, k].velocity.x;
							newVelocity.y = -1 *cells[i, j-1, k].velocity.y;
							newVelocity.z = 	cells[i, j-1, k].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i, j-1, k].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (front) {
							newVelocity.x = 	cells[i, j, k+1].velocity.x;
							newVelocity.y = 	cells[i, j, k+1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k+1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i, j, k+1].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						} else if (back) {
							newVelocity.x = 	cells[i, j, k-1].velocity.x;
							newVelocity.y =		cells[i, j, k-1].velocity.y;
							newVelocity.z = -1 *cells[i, j, k-1].velocity.z;
							cells[i, j, k].SetVelocity(newVelocity);
							//cells[i, j, k].SetDensity(cells[i, j, k-1].density);
							cells[i, j, k].SetDensity(0f, 0);
							cells[i, j, k].SetDensity(0f, 1);
							cells[i, j, k].SetDensity(0f, 2);
						}
						Debug.DrawLine(cells[i, j, k].position, cells[i, j, k].position + cells[i, j, k].velocity, Color.blue);
					}
				}
			}
		}
	}
}