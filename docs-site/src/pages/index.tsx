import Link from '@docusaurus/Link';
import Layout from '@theme/Layout';
import useBaseUrl from '@docusaurus/useBaseUrl';
import styles from './index.module.css';

export default function Home() {
  const logoUrl = useBaseUrl('/img/cbicall-logo-transparent.png');

  return (
    <Layout
      title="CBIcall"
      description="CNAG Biomedical Informatics framework for variant calling">
      <main>
        <section className={styles.hero}>
          <div className={styles.heroGrid}>
            <div className={styles.heroCopy}>
              <img className={styles.heroLogo} src={logoUrl} alt="CBIcall" />
              <h1 className={styles.srOnly}>CBIcall</h1>
              <p className={styles.heroSubtitle}>
                CBIcall is a configuration-driven framework for reproducible
                variant calling in large sequencing cohorts.
              </p>
              <div className={styles.actions}>
                <Link className="button button--primary button--lg" to="/docs/usage/quickstart">
                  Quick Start
                </Link>
                <Link className="button button--secondary button--lg" to="/docs/overview">
                  Overview
                </Link>
              </div>
            </div>

            <div className={styles.workflowPanel} aria-label="CBIcall framework summary">
              <div className={styles.workflowStep}>
                <span>Parameters</span>
                <strong>Validated YAML analysis contract</strong>
              </div>
              <div className={styles.workflowArrow}>→</div>
              <div className={styles.workflowStepPrimary}>
                <span>Workflow backends</span>
                <strong>Bash · Snakemake · Nextflow · Cromwell</strong>
              </div>
              <div className={styles.workflowArrow}>→</div>
              <div className={styles.workflowStep}>
                <span>Audit</span>
                <strong>Run reports · hashes · comparisons</strong>
              </div>
              <div className={styles.tokens}>
                <span>validation</span>
                <span>resources</span>
                <span>provenance</span>
                <span>reports</span>
                <span>nf-core</span>
              </div>
            </div>
          </div>
        </section>
        <section className={styles.sections}>
          <div className={styles.grid}>
            <Link className={styles.card} to="/docs/installation/docker">
              <span>Setup</span>
              <h2>Install</h2>
              <p>Docker, Apptainer, and non-containerized setup paths.</p>
            </Link>
            <Link className={styles.card} to="/docs/pipelines/overview">
              <span>Pipelines</span>
              <h2>Included Analyses</h2>
              <p>Curated WES, WGS, and mtDNA pipelines for single-sample and cohort use.</p>
            </Link>
            <Link className={styles.card} to="/docs/help/run-comparison">
              <span>Reproducibility</span>
              <h2>Audit Runs</h2>
              <p>Structured reports, workflow fingerprints, resource hashes, and run comparison.</p>
            </Link>
            <Link className={styles.card} to="/docs/technical-details/adding-a-pipeline">
              <span>Developers</span>
              <h2>Extend</h2>
              <p>Registry-driven backend wiring for adding new pipelines.</p>
            </Link>
          </div>
        </section>
      </main>
    </Layout>
  );
}
